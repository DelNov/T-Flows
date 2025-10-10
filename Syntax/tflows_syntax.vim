"-------------------------------------
" Tuning the color scheme for T-Flows
"-------------------------------------

" Keyords should be Fortran commands
highlight fortranKeyword ctermfg=Yellow

" Ghost numbers are not desirable, mark them red
highlight fortranNumber  ctermfg=Red
highlight fortranFloat  ctermfg=Red

" Constants are better than ghost numbers, mark them blue ...
" ... they are safer, they are in effect former ghost numbers
highlight fortranConstant  ctermfg=Blue

" Continuation lines.  I like them visible and a bit psychedelic
highlight fortranContinueMark  ctermfg=Magenta

" Strings ... they are kind of lame, paint them gray
highlight fortranString  ctermfg=Gray

" Operators ... I don't know, maybe the same as keywords
highlight fortranOperator  ctermfg=Yellow

" Just Comment didn't work, fortranComment did
highlight fortranComment  ctermfg=LightBlue

" Try to over-ride Fortran objects in T-Flows
highlight fortranObjectTflows cterm=Bold ctermfg=Green

" MPI calls in T-Flows - they are few and far in between, yet dangerous
highlight fortranMpiTflows  ctermfg=Red

" OpenACC links in T-Flows
highlight fortranGPUTflows  ctermfg=LightMagenta

" Calls to PETSc functions which are written in C
" I am not even sure they do the argument checking
highlight fortranPetscTflows  ctermfg=LightYellow

" Macros in T-Flows, maybe the same color as Fortran keywords, but bold
highlight fortranMacroTflows  cterm=Bold ctermfg=Yellow

