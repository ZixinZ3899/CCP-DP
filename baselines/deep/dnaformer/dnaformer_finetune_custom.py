import argparse
import csv
import json
import random
import time
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import torch
import torch.nn.functional as F
from torch.utils.data import DataLoader, Dataset


ADAPTER_DIR = Path(__file__).resolve().parent
REPO_ROOT = ADAPTER_DIR.parents[2]
DNAFORMER_ROOT = (
    REPO_ROOT
    / "build"
    / "third_party"
    / "TReconLM"
    / "DeepLearningBaselines"
    / "DNAFormer"
)
sys.path.insert(0, str(DNAFORMER_ROOT))

from model_DNAFormer_siamese import net


DNA_TO_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}
INDEX_TO_DNA = "ACGT"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--train", required=True)
    parser.add_argument("--val", required=True)
    parser.add_argument("--test", required=True)
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--target_len", type=int, required=True)
    parser.add_argument("--input_len", type=int, default=132)
    parser.add_argument("--max_reads", type=int, default=10)
    parser.add_argument("--epochs", type=int, default=30)
    parser.add_argument("--batch_size", type=int, default=8)
    parser.add_argument("--accumulation_steps", type=int, default=8)
    parser.add_argument("--backbone_lr", type=float, default=1e-5)
    parser.add_argument("--new_layer_lr", type=float, default=1e-3)
    parser.add_argument("--weight_decay", type=float, default=0.0)
    parser.add_argument("--patience", type=int, default=6)
    parser.add_argument("--seed", type=int, default=2026)
    parser.add_argument("--num_workers", type=int, default=0)
    return parser.parse_args()


def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def encode_indices(sequence):
    return torch.tensor([DNA_TO_INDEX[x] for x in sequence], dtype=torch.long)


def encode_read(sequence):
    return F.one_hot(encode_indices(sequence), num_classes=4).transpose(0, 1).float()


class ClusterDataset(Dataset):
    def __init__(self, path, target_len, input_len, max_reads, training, seed):
        self.path = Path(path)
        self.target_len = target_len
        self.input_len = input_len
        self.max_reads = max_reads
        self.training = training
        self.seed = seed
        self.examples = []

        with self.path.open() as handle:
            for line_number, raw in enumerate(handle, 1):
                line = raw.strip()
                if not line:
                    continue
                reads_text, center = line.rsplit(":", 1)
                reads = reads_text.split("|")
                if len(center) != target_len:
                    raise ValueError(
                        f"{self.path}:{line_number}: center length "
                        f"{len(center)} != {target_len}"
                    )
                if len(reads) < 2:
                    raise ValueError(f"{self.path}:{line_number}: fewer than 2 reads")
                if set(center) - set(DNA_TO_INDEX):
                    raise ValueError(f"{self.path}:{line_number}: invalid center")
                for read in reads:
                    if set(read) - set(DNA_TO_INDEX):
                        raise ValueError(f"{self.path}:{line_number}: invalid read")
                    if len(read) > input_len:
                        raise ValueError(
                            f"{self.path}:{line_number}: read length "
                            f"{len(read)} > {input_len}"
                        )
                self.examples.append((reads, center))

    def __len__(self):
        return len(self.examples)

    def _select_reads(self, reads, index):
        count = min(self.max_reads, len(reads))
        if len(reads) <= self.max_reads:
            selected = list(reads)
            if self.training:
                random.shuffle(selected)
            return selected
        if self.training:
            return random.sample(reads, count)
        rng = random.Random(self.seed + 1_000_003 * index)
        chosen = rng.sample(range(len(reads)), count)
        return [reads[i] for i in chosen]

    def __getitem__(self, index):
        reads, center = self.examples[index]
        selected = self._select_reads(reads, index)

        left = torch.zeros(self.max_reads, 4, self.input_len, dtype=torch.float32)
        right = torch.zeros_like(left)

        for read_index, read in enumerate(selected):
            encoded = encode_read(read)
            length = encoded.shape[-1]
            left[read_index, :, :length] = encoded
            right[read_index, :, :length] = torch.flip(encoded, dims=[-1])

        return left, right, encode_indices(center)


def build_config(args):
    return SimpleNamespace(
        label_length=args.target_len,
        filter_index=False,
        index_length=0,
        noisy_copies_length=args.input_len,
        max_number_per_cluster=args.max_reads,
        n_head=32,
        activation="gelu",
        num_layers=12,
        d_model=1024,
        alignment_filters=128,
        dim_feedforward=2048,
        output_ch=4,
        enc_filters=4,
        p_dropout=0.0,
        class_token=0,
        use_input_scaling=False,
        device="cuda",
    )


def clean_official_state(checkpoint):
    raw = checkpoint.get("model_state_dict", checkpoint.get("state_dict", checkpoint))
    cleaned = {}
    for old_key, value in raw.items():
        key = old_key
        for prefix in ("module._orig_mod.", "module.", "_orig_mod."):
            if key.startswith(prefix):
                key = key[len(prefix):]
                break
        cleaned[key] = value
    if "fusion.pred_fusion" in cleaned:
        cleaned["fusion.pred_fusion_left"] = cleaned.pop("fusion.pred_fusion")
    if "fusion.pred_fusion_flip" in cleaned:
        cleaned["fusion.pred_fusion_right"] = cleaned.pop("fusion.pred_fusion_flip")
    return cleaned


def load_transfer_weights(model, checkpoint_path):
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=False)
    official = clean_official_state(checkpoint)
    current = model.state_dict()
    compatible = {
        key: value
        for key, value in official.items()
        if key in current and tuple(value.shape) == tuple(current[key].shape)
    }
    result = model.load_state_dict(compatible, strict=False)
    if result.unexpected_keys:
        raise RuntimeError(f"Unexpected checkpoint keys: {result.unexpected_keys}")
    new_names = set(result.missing_keys)
    print("Checkpoint epoch:", checkpoint.get("epoch", "unknown"))
    print("Transferred tensors:", len(compatible))
    print("New tensors:", len(new_names))
    for name in sorted(new_names):
        print("  new:", name, tuple(current[name].shape))
    return new_names


def calculate_loss(outputs, labels):
    label_one_hot = F.one_hot(labels, num_classes=4).permute(0, 2, 1).float()
    main = F.cross_entropy(outputs["pred"], labels)
    left = F.cross_entropy(outputs["pred_left"].softmax(dim=1), label_one_hot)
    right = F.cross_entropy(outputs["pred_right"].softmax(dim=1), label_one_hot)
    return main + 0.25 * 0.5 * (left + right)


def edit_distance(a, b):
    if len(a) > len(b):
        a, b = b, a
    previous = list(range(len(a) + 1))
    for i, char_b in enumerate(b, 1):
        current = [i]
        for j, char_a in enumerate(a, 1):
            current.append(
                min(
                    current[-1] + 1,
                    previous[j] + 1,
                    previous[j - 1] + (char_a != char_b),
                )
            )
        previous = current
    return previous[-1]


def indices_to_sequence(indices):
    return "".join(INDEX_TO_DNA[x] for x in indices)


@torch.inference_mode()
def evaluate(model, loader, device, timed=False):
    model.eval()
    if timed:
        torch.cuda.synchronize()
    start = time.perf_counter()

    total_loss = 0.0
    total_examples = 0
    correct_bases = 0
    total_bases = 0
    predictions = []
    references = []

    for left, right, labels in loader:
        left = left.to(device, non_blocking=True)
        right = right.to(device, non_blocking=True)
        labels = labels.to(device, non_blocking=True)
        model_input = torch.cat([left, right], dim=0)

        with torch.cuda.amp.autocast(dtype=torch.float16):
            outputs = model(model_input)
            loss = calculate_loss(outputs, labels)

        predicted = outputs["pred"].argmax(dim=1)
        correct_bases += predicted.eq(labels).sum().item()
        total_bases += labels.numel()
        total_examples += labels.shape[0]
        total_loss += loss.item() * labels.shape[0]
        predictions.extend(predicted.cpu().tolist())
        references.extend(labels.cpu().tolist())

    if timed:
        torch.cuda.synchronize()
    seconds = time.perf_counter() - start

    predicted_strings = [indices_to_sequence(x) for x in predictions]
    reference_strings = [indices_to_sequence(x) for x in references]
    distances = [edit_distance(p, r) for p, r in zip(predicted_strings, reference_strings)]
    success_count = sum(distance == 0 for distance in distances)
    total_ed = sum(distances)

    return {
        "loss": total_loss / total_examples,
        "clusters": total_examples,
        "base_accuracy": correct_bases / total_bases,
        "success_count": success_count,
        "success_rate": success_count / total_examples,
        "total_edit_distance": total_ed,
        "average_edit_distance": total_ed / total_examples,
        "seconds": seconds,
        "clusters_per_second": total_examples / seconds,
        "predictions": predicted_strings,
    }


def metrics_without_predictions(metrics):
    return {key: value for key, value in metrics.items() if key != "predictions"}


def main():
    args = parse_args()
    set_seed(args.seed)

    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is required")
    device = torch.device("cuda")
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    train_dataset = ClusterDataset(
        args.train, args.target_len, args.input_len, args.max_reads, True, args.seed
    )
    val_dataset = ClusterDataset(
        args.val, args.target_len, args.input_len, args.max_reads, False, args.seed + 10_000
    )
    test_dataset = ClusterDataset(
        args.test, args.target_len, args.input_len, args.max_reads, False, args.seed + 20_000
    )

    generator = torch.Generator().manual_seed(args.seed)
    loader_kwargs = dict(
        batch_size=args.batch_size,
        num_workers=args.num_workers,
        pin_memory=True,
    )
    train_loader = DataLoader(
        train_dataset, shuffle=True, generator=generator, **loader_kwargs
    )
    val_loader = DataLoader(val_dataset, shuffle=False, **loader_kwargs)
    test_loader = DataLoader(test_dataset, shuffle=False, **loader_kwargs)

    print("Device:", device)
    print("GPU:", torch.cuda.get_device_name(0))
    print("Train/val/test:", len(train_dataset), len(val_dataset), len(test_dataset))
    print("Batch size:", args.batch_size)
    print("Accumulation steps:", args.accumulation_steps)
    print("Effective batch size:", args.batch_size * args.accumulation_steps)

    model = net(build_config(args))
    new_names = load_transfer_weights(model, args.checkpoint)
    model = model.to(device)

    backbone_parameters = []
    new_parameters = []
    for name, parameter in model.named_parameters():
        if name in new_names:
            new_parameters.append(parameter)
        else:
            backbone_parameters.append(parameter)

    optimizer = torch.optim.AdamW(
        [
            {"params": backbone_parameters, "lr": args.backbone_lr},
            {"params": new_parameters, "lr": args.new_layer_lr},
        ],
        weight_decay=args.weight_decay,
    )
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=2, min_lr=1e-7
    )
    scaler = torch.cuda.amp.GradScaler()

    best_path = output_dir / "best_model.pth"
    log_path = output_dir / "training_log.csv"
    best_key = None
    best_epoch = 0
    epochs_without_improvement = 0

    with log_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "epoch",
                "train_loss",
                "train_base_accuracy",
                "val_loss",
                "val_base_accuracy",
                "val_success_count",
                "val_success_rate",
                "val_total_edit_distance",
                "val_average_edit_distance",
                "backbone_lr",
                "new_layer_lr",
                "seconds",
            ]
        )

        torch.cuda.reset_peak_memory_stats()
        training_start = time.perf_counter()

        for epoch in range(1, args.epochs + 1):
            model.train()
            epoch_start = time.perf_counter()
            optimizer.zero_grad(set_to_none=True)
            train_loss_sum = 0.0
            train_examples = 0
            train_correct_bases = 0
            train_total_bases = 0

            for batch_index, (left, right, labels) in enumerate(train_loader):
                left = left.to(device, non_blocking=True)
                right = right.to(device, non_blocking=True)
                labels = labels.to(device, non_blocking=True)
                model_input = torch.cat([left, right], dim=0)

                with torch.cuda.amp.autocast(dtype=torch.float16):
                    outputs = model(model_input)
                    raw_loss = calculate_loss(outputs, labels)
                    loss = raw_loss / args.accumulation_steps

                scaler.scale(loss).backward()

                should_step = (
                    (batch_index + 1) % args.accumulation_steps == 0
                    or (batch_index + 1) == len(train_loader)
                )
                if should_step:
                    scaler.unscale_(optimizer)
                    torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
                    scaler.step(optimizer)
                    scaler.update()
                    optimizer.zero_grad(set_to_none=True)

                predicted = outputs["pred"].detach().argmax(dim=1)
                train_correct_bases += predicted.eq(labels).sum().item()
                train_total_bases += labels.numel()
                train_examples += labels.shape[0]
                train_loss_sum += raw_loss.item() * labels.shape[0]

            val_metrics = evaluate(model, val_loader, device)
            scheduler.step(val_metrics["average_edit_distance"])
            epoch_seconds = time.perf_counter() - epoch_start
            train_loss = train_loss_sum / train_examples
            train_base_accuracy = train_correct_bases / train_total_bases
            current_key = (
                val_metrics["total_edit_distance"],
                -val_metrics["success_count"],
                -val_metrics["base_accuracy"],
            )

            if best_key is None or current_key < best_key:
                best_key = current_key
                best_epoch = epoch
                epochs_without_improvement = 0
                torch.save(
                    {
                        "epoch": epoch,
                        "model_state_dict": model.state_dict(),
                        "val_metrics": metrics_without_predictions(val_metrics),
                        "args": vars(args),
                        "transferred_from": str(args.checkpoint),
                    },
                    best_path,
                )
                marker = "BEST"
            else:
                epochs_without_improvement += 1
                marker = ""

            row = [
                epoch,
                train_loss,
                train_base_accuracy,
                val_metrics["loss"],
                val_metrics["base_accuracy"],
                val_metrics["success_count"],
                val_metrics["success_rate"],
                val_metrics["total_edit_distance"],
                val_metrics["average_edit_distance"],
                optimizer.param_groups[0]["lr"],
                optimizer.param_groups[1]["lr"],
                epoch_seconds,
            ]
            writer.writerow(row)
            handle.flush()

            print(
                f"epoch={epoch}/{args.epochs} "
                f"train_loss={train_loss:.6f} "
                f"train_base={train_base_accuracy:.6f} "
                f"val_base={val_metrics['base_accuracy']:.6f} "
                f"val_success={val_metrics['success_count']}/{val_metrics['clusters']} "
                f"({val_metrics['success_rate']:.6f}) "
                f"val_avg_ed={val_metrics['average_edit_distance']:.6f} "
                f"seconds={epoch_seconds:.2f} {marker}"
            )

            if epochs_without_improvement >= args.patience:
                print("Early stopping at epoch", epoch)
                break

        training_seconds = time.perf_counter() - training_start
        peak_training_gpu_mib = torch.cuda.max_memory_allocated() / 1024**2

    best_checkpoint = torch.load(best_path, map_location="cpu", weights_only=False)
    model.load_state_dict(best_checkpoint["model_state_dict"], strict=True)
    model.to(device)

    torch.cuda.reset_peak_memory_stats()
    test_metrics = evaluate(model, test_loader, device, timed=True)
    peak_inference_gpu_mib = torch.cuda.max_memory_allocated() / 1024**2

    predictions_path = output_dir / "test_predictions.txt"
    with predictions_path.open("w") as handle:
        for sequence in test_metrics["predictions"]:
            handle.write(sequence + "\n")

    final_metrics = metrics_without_predictions(test_metrics)
    final_metrics.update(
        {
            "best_epoch": best_epoch,
            "training_seconds": training_seconds,
            "peak_training_gpu_memory_mib": peak_training_gpu_mib,
            "peak_inference_gpu_memory_mib": peak_inference_gpu_mib,
            "best_checkpoint": str(best_path.resolve()),
            "predictions_file": str(predictions_path.resolve()),
            "test_file": str(Path(args.test).resolve()),
        }
    )

    metrics_path = output_dir / "test_metrics.json"
    with metrics_path.open("w") as handle:
        json.dump(final_metrics, handle, indent=2)

    print("=" * 72)
    print("Best epoch:", best_epoch)
    print("Training seconds:", training_seconds)
    print("Peak training GPU memory MiB:", peak_training_gpu_mib)
    print(json.dumps(final_metrics, indent=2))
    print("DNAFORMER FINE-TUNING: SUCCESS")


if __name__ == "__main__":
    main()