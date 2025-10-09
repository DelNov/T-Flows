"----------------------------------
" Tuning the color scheme for Gmsh
"----------------------------------

" (1) If, For, Function are all the same color, bold green
highlight gmshConditional cterm=Bold ctermfg=Green
highlight gmshRepeat      cterm=Bold ctermfg=Green
highlight gmshFunction    cterm=Bold ctermfg=Green

" (2) String
highlight gmshString ctermfg=Gray

" (3) Numbers (ghosts) are not desirable, mark them red
highlight gmshNumber ctermfg=Red
highlight gmshFloat  ctermfg=Red

" (4) Constants
highlight gmshConstant ctermfg=Blue

" (5) Comment
highlight gmshComment ctermfg=Magenta

" (6) Keyword
highlight gmshKeyword ctermfg=Yellow
