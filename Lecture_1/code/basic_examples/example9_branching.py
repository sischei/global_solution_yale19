def hat_function(x):
  """hat function
  """
  if x < 0:
    return 0.0
  elif 0 <= x < 1:
    return x
  elif 1 <= x < 2:
    return 2 - x
  elif x >= 2:
    return 0.0
  
print hat_function(1.4)
