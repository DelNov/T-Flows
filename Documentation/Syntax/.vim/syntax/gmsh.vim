" Vim syntax file
" Language:	gmsh
"

" Remove any old syntax stuff hanging around
set tabstop=2
set shiftwidth=2
set expandtab
syn clear

"-----------------------
" (1) If, For, Function
"-----------------------
syn keyword gmshConditional If Else EndIf
syn keyword gmshRepeat      For EndFor
syn keyword gmshFunction    Function Macro Return

"-------------
" (2) Strings
"-------------
syn region gmshString       start=+'+ end=+'+ oneline
syn region gmshString       start=+"+ end=+"+ oneline

"--------------------------------
" (3) Numbers of different sorts
"--------------------------------

" Integer
syn match gmshNumber        "-\=\<\d\+\>"

" Floating point number, with dot, optional exponent
syn match gmshFloat         "\<\d\+\.\d*\([edED][-+]\=\d\+\)\=\>"
" Floating point number, starting with a dot, optional exponent
syn match gmshFloat         "\.\d\+\([edED][-+]\=\d\+\)\=\>"
" Floating point number, without dot, with exponent
syn match gmshFloat         "\<\d\+[edED][-+]\=\d\+\>"

"---------------
" (4) Constants
"---------------
syn match gmshConstant      "\<[A-Z_][A-Z0-9_]*\>"

"--------------
" (5) Comments
"--------------
syn region gmshComment      start="//" skip="\\$" end="$" keepend contains=@gmshCommentGroup,gmshSpaceError
syn region gmshComment      matchgroup=gmshCommentStart start="/\*" matchgroup=NONE end="\*/" contains=@gmshCommentGroup,gmshCommentStartError,gmshSpaceError

"--------------
" (6) Keywords
"--------------
syn keyword gmshKeyword newreg newp Abort Acos ArcCos Asin ArcSin Atan ArcTan
syn keyword gmshKeyword Atan2 ArcTan2 Attractor Bump BSpline Bounds BoundingBox
syn keyword gmshKeyword Ceil Cosh Cos Characteristic Circle Coherence Complex
syn keyword gmshKeyword Color ColorTable CatmullRom Call Curve Delete Dilate
syn keyword gmshKeyword Duplicata Draw Exp Ellipsis Extrude Elliptic
syn keyword gmshKeyword Exit Fabs Floor Fmod
syn keyword gmshKeyword Hypot In Include Intersect Knots Length
syn keyword gmshKeyword Line Loop Log Log10 Layers Modulo Meshes Nurbs
syn keyword gmshKeyword Order Physical Pi Plane Point Power Progression
syn keyword gmshKeyword Parametric Printf Recombine Rotate Ruled Rand
syn keyword gmshKeyword Sqrt Sin Sinh Spline Surface Symmetry
syn keyword gmshKeyword Sprintf Transfinite Translate Tanh Tan Trimmed
syn keyword gmshKeyword Using Volume With SS VS TS ST VT TT SL VL TL SP VP TP
syn keyword gmshKeyword System

if !exists("did_gmsh_syntax_inits")
  let did_gmsh_syntax_inits = 1

  " The default methods for highlighting.  Can be overridden later

  " (1) If, For, Function
  hi link gmshConditional Conditional
  hi link gmshRepeat      Repeat
  hi link gmshFunction    Function

  " (2) String
  hi link gmshString      String

  " (3) All sorts of numbers
  hi link gmshNumber      Number
  hi link gmshFloat       Float

  " (4) Constants
  hi link gmshConstant    Constant

  " (5) Comments
  hi link gmshComment     Comment

  " (6) Keyword
  hi link gmshKeyword     Keyword
endif

let b:current_syntax = "gmsh"

