# Tophat port to Agave

Notes
------
* Issue: Tophat app relies on presence of iget to grab all the child index files associated with a given reference
    * Addressed. Now uses refprep-based bundles
* Tophat bundles its own binaries. This is good, but there are a handful of module dependencies. We want to get rid of those to make the app portable to VMs in the Cloud
* Issue: Relies on system samtools, and is of an older vintage so need to pin to samtool-0.1.9 instead of 1.x series
    * Addressed
* Implements a nice cleanup section.
    * This needed to be updated to accomodate the tar format for reference genomes
* Issue: Uses explicit parameters in the bash script. Would be good to port this to showArguments logic in the future to simplify things

Log
```
auth-tokens-create -S
API password:
Token for iplantc.org:vaughn successfully refreshed and cached for 14400 seconds
f651c4e2a391d014226dadc279096f4

files-upload -F dnalc-tophat-hpc-2.0.11.1 vaughn/applications

apps-addupdate -F dnalc-tophat-hpc-2.0.11.1/tophat_lonestar.json
Successfully added app dnalc-tophat-lonestar-2.0.11.1

jobs-submit -F dnalc-tophat-lonestar-job.json -W
jobs-history -v 0001428683795837-5056a550b8-0001-007

```