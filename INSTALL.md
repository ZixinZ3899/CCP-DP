# Install into the CCP-DP repository

Extract this archive from the root of the existing CCP-DP repository. It adds
documentation and scripts under `baselines/` and does not contain datasets,
results, compiled programs, or third-party source already present in the
repository.

```bash
cd ~/CCP-DP
unzip -o /path/to/CCPDP_baseline_readmes_and_scripts.zip
chmod +x baselines/build_baselines.sh baselines/itr/fetch_upstream.sh
```

Then run:

```bash
./baselines/itr/fetch_upstream.sh
./baselines/build_baselines.sh
```

Read `baselines/README.md` for complete dataset and execution examples.

