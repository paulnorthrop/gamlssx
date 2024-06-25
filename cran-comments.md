## Resubmission

In this version I have 

* DESCRIPTION: removed unnecessary spaces at the ends of lines
* fitGEV.R: avoided writing messages to the console. By default no messages are written. The user can turn them on.
* fitGEV.R: in ... noted that trace = FALSE can be used to avoid gamlss() writing to the console. 

## R CMD check results

0 errors | 0 warnings | 0 notes

The possibly misspelled words in DESCRIPTION: 
  Rigby (14:14), 
  Stasinopoulos (14:30)
are false positives (author names)
  
## Test environments

- macOS (R-release), ubuntu (R-oldrel, R-release, R-devel), windows (R-release) using the rcmdcheck package
- macOS builder 
- win-builder (R-devel, R-release and R-oldrelease)

## Downstream dependencies

None, this is a new submission
