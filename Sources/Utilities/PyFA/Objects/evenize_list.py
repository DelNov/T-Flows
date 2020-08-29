
def evenize_list(low, high, n):

  numbers = []

  for i in range(1, n+1):
    numbers.append( int(low + i  * (high-low) / (n+1)) )

  return numbers

