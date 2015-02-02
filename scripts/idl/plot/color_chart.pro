index = [1, 3, 4, 7, 9, 13, 14, 16, 18, 19, 20, 21, 22, 24, 25, 26]
FOR j=0,15 DO cgLoadCT, index[j], NCOLORS=16, BOTTOM=16*j, /BREWER
CIndex, /BREWER
