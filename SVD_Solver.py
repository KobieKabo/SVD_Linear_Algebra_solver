import numpy as np

def svd_solver(A):
  # Grabs m & n size from matrix A
  m, n = A.shape
  v = np.random.rand(n)
  u = np.dot(A.T, v)
  u_norm = np.linalg.norm(u)
  u = u / u_norm

  v = np.dot(A, u)
  v_norm = np.linalg.norm(v)
  v = v / v_norm

  sigma = np.linalg.norm(np.dot(A, v))
  u = np.dot(A.T, v) / sigma

  return u, sigma, v

def main():
  rows = int(input("Enter the number of rows in matrix A: "))
  cols = int(input("Enter the number of columns in matrix A: "))
  A = np.zeros((rows, cols))

  for i in range(rows):
    for j in range(cols):
      A[i, j] = float(input(f"Enter element A[{i+1},{j+1}]: "))

  u, sigma, v = svd_solver(A)

  print("\nU:")
  print(u)
  print("\nSigma:")
  print(sigma)
  print("\nV:")
  print(v)

if __name__=="__main__":
  main()