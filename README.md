# dnalc_fAPI_apps
DNALC's iPlant foundation API apps

* dnalc-fastqc-stampede-0.10.1u1  @jamescarson3
* dnalc-fxtrim-stampede-0.0.13.3u1	@jamescarson3
* dnalc-tophat-stampede-2.0.11.1u1    @mwvaughn
* dnalc-cuffdiff-stampede-2.1.1.2u1   @johnfonner
* dnalc-cufflinks-stampede-2.1.1.2u1  @johnfonner
* dnalc-cuffmerge-stampede-2.1.1u1    @johnfonner

Porting guidelines
* Deploy to a Lonestar4 clone. Use normal queue
* Use showArguments where possible
* Make sure to preserve the original logic and interfaces of the original wrapper scripts
* Create a test job for each application, including input and output files
