# Tophat port to Agave

Notes
------
* Issue: Tophat app relies on presence of iget to grab all the child index files associated with a given reference
* Tophat bundles its own binaries. This is good
* Issue: Relies on system samtools, and is of an older vintage so need to pin to samtool-0.1.9 instead of 1.x series
* Implements a nice cleanup section
* Issue: Uses explicit parameters in the bash script. Would be good to port this to showArguments logic in the future to simplify things

Log
```
auth-tokens-create -S
API password:
Token for iplantc.org:vaughn successfully refreshed and cached for 14400 seconds
f651c4e2a391d014226dadc279096f4

files-upload -F dnalc-tophat-hpc-2.0.8 vaughn/applications

apps-addupdate -F dnalc-tophat-hpc-2.0.8/tophat_lonestar.json

```