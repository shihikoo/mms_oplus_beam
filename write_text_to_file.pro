
FUNCTION write_text_to_file, filename, text, APPEND = APPEND
  OPENU, unit, filename, APPEND = APPEND,/GET_LUN 
  PRINTF, unit, text
  FREE_LUN, unit   
  close, /all
END
