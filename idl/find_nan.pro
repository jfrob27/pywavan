;This function find NaNs in an array and replace them with zeros

FUNCTION find_nan, array

nan = where(finite(array) EQ 0)

array[nan] = 0.

return, array

END