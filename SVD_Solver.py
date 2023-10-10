import numpy as np
# Jakob Long, JRL4725

def svd_solver(A):
    """Calculate SVD for a general NXM matrix.
    
    Args:
        A : Input matrix of shape (m, n).
        
    Returns:
        U : Left singular vectors of shape (m, m).
        S : Singular values in descending order of shape (min(m, n),).
        V : Right singular vectors of shape (n, n).
        condition_number : Condition number of the matrix.
        A_inv : Inverse of the matrix A.
    """
    m,n = A.shape
    
    # Compute eigenvalues and eigenvectors of U and V
    eigen_val_U, eigen_vec_U = np.linalg.eigh(np.dot(A,A.T))
    eigen_val_V, eigen_vec_V = np.linalg.eigh(np.dot(A.T,A))
    
    # Calculate condition number
    sigma = np.sqrt(eigen_val_U)
    condition_number = (max(sigma) / min(sigma))
    # Compute Sigma
    sigma = np.zeros((m,n))
    for i in range(min(len(eigen_val_U),len(eigen_val_V))):
      sigma[i,i] = np.sqrt(eigen_val_U[i])
    
    # Sort eigenval & eigenvec 
    sort_U = np.argsort(eigen_val_U)[::-1]
    sort_V = np.argsort(eigen_val_V)[::-1] 
    
    eigen_val_U = eigen_val_U[sort_U]
    U = eigen_vec_U[:,sort_U]
    
    eigen_val_V = eigen_val_V[sort_V]
    V = eigen_vec_V[:,sort_V]
    
    # Calculate inverse
    try:
      for i in range(min(n,m)):
        if sigma[i][i] == 0:
          raise Exception(
            'Error: Matrix is singular & has no inverse.')
      if (np.diag(sigma)).any() == 0:
        raise Exception(
          'Error: Matrix is singular & has no inverse.')
      else:
        sigma_inv = np.zeros(n,m)
        for i in range(min(n,m)):
          sigma_inv[i][i] == 1/sigma[i][i]
        A_inv = V @ sigma_inv @ U.T
    except:
      A_inv = 'Error: Matrix is singular & has no inverse.'
    
    return U, sigma, V, condition_number, A_inv

def main():
  while True:
    try:
      rows = int(input("Enter the number of rows in matrix A: "))
      cols = int(input("Enter the number of columns in matrix A: "))
      if rows <= 0 or cols <= 0:
        print('Please input a valid integer, of 1 or larger.')
        continue
      break
    except ValueError:
      print('Enter a valid integer of 1 or larger.')
      continue
  
  A = np.zeros((rows, cols))

  for i in range(rows):
    for j in range(cols):
      A[i, j] = float(input(f"Enter element A[{i+1},{j+1}]: "))

  u, s, v, cond_num, A_inv = svd_solver(A)
  u1,s1,v1 = np.linalg.svd(A)
  
  print("\nA:")
  print(A)
  print('\n___________________________')
  print('\nA Reconstruction:')
  print(u@s@v.T)
  print('\n___________________________')
  print("\nU:")
  print(u)
  print('\n___________________________')
  print("\nSigma:")
  print(s)
  print('\n___________________________')
  print("\nV:")
  print(v.T)
  print('\n___________________________')
  print("\nCondition Numbers:")
  print(cond_num)
  print('\n___________________________')
  print("\nA-Inverse:")
  print(A_inv)
  print('\n___________________________')
  print("\nBlackbox U:")
  print(u1)
  print('\n___________________________')
  print("\nBlackBox Sigma:")
  print(s1)
  print('\n___________________________')
  print("\nBlackBox V:")
  print(v1)
  print('\n___________________________')
  

if __name__=="__main__":
  main()