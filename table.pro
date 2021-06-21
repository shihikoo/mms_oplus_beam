;-----------------------------------------------
; Purpose: The function export the unique values in a vector and
; number of occurance of that value.
; Inputs: vec, a vector
;
; Written by Jing Liao
; Written on 05/14/2021
;-----------------------------------------------

FUNCTION table, vec
  values = vec[uniq(vec,sort(vec))]
  nvalues = N_ELEMENTS(values)
  output = STRARR(nvalues, 2)
  FOR i = 0, nvalues-1 DO BEGIN 
     value = values(i)
     index = where(vec EQ value, ct)
     output[i, 0] = value
     output[i, 1] = ct
  ENDFOR
  RETURN,output
END
