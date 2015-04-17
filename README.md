# dnalc_fAPI_apps
DNALC's iPlant foundation API apps

# Original app names
------------------
* dnalc-fastqc-stampede-0.10.1u1  @jamescarson3
* dnalc-fxtrim-stampede-0.0.13.3u1	@jamescarson3
* dnalc-tophat-stampede-2.0.11.1u1    @mwvaughn
* dnalc-cuffdiff-stampede-2.1.1.2u1   @johnfonner
* dnalc-cufflinks-stampede-2.1.1.2u1  @johnfonner
* dnalc-cuffmerge-stampede-2.1.1u1    @johnfonner

# New app names
------------------
* dnasubway-fastqc-lonestar-0.11.2.0
* dnasubway-fastx-lonestar-0.0.13.2.0
* dnasubway-tophat-refprep-lonestar-2.0.11.1
* dnasubway-tophat-lonestar-2.0.11.1
* dnasubway-cufflinks-lonestar-2.1.1
* dnasubway-cuffdiff-lonestar-2.1.1
* dnasubway-cuffmerge-lonestar-2.1.1

Porting guidelines
------------------
* Deploy to a Lonestar4 clone. Use normal queue
* Use showArguments where possible
* Make sure to preserve the original logic and interfaces of the original wrapper scripts
* Create a test job for each application, including input and output files
    * /iplant/home/shared/iplant_DNA_subway is writeable by all project maintainers
