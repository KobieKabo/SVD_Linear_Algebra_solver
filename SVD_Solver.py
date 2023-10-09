import numpy as np

def svd_solver(A):
    """Calculate SVD for a general NXM matrix.
    
    Args:
        A : Input matrix of shape (m, n).
        
    Returns:
        U : Left singular vectors of shape (m, m).
        S : Singular values in descending order of shape (min(m, n),).
        Vt : Right singular vectors of shape (n, n).
        condition_number (float): Condition number of the matrix.
        A_inv : Inverse of the matrix A.
    """
    m, n = A.shape
    
    # Step 2: Compute eigenvalues and eigenvectors of ATA and AAT
    eigen_val_U, eigen_vec_U = np.linalg.eigh(np.dot(A,A.T))
    eigen_val_V, eigen_vec_V = np.linalg.eigh(np.dot(A.T,A))
    
    # Step 3: Sort eigenvalues in descending order
    eigvals_U = np.flip(eigen_val_U)
    eigvecs_U = np.fliplr(eigen_vec_U)
    eigvals_V = np.flip(eigen_val_V)
    eigvecs_V = np.fliplr((eigen_vec_V))
    
    # Step 4: Compute singular values and sort by magnitude
    sigma = np.sqrt(eigvals_V)
    singular_value_indices = np.argsort(sigma)[::-1]
    sigma = sigma[singular_value_indices]
    
    # Step 5: Compute U, S, and Vt
    U = eigvecs_U[:, singular_value_indices]
    Vt = eigvecs_V.T[:, singular_value_indices]
    S = np.zeros((m, n))
    np.fill_diagonal(S, sigma)
    
    # Step 6: Calculate condition number
    condition_number = np.max(sigma) / np.min(sigma)
    
    # Step 7: Calculate pseudo-inverse
    S_inv = np.linalg.inv(S)
    try:
      A_inv = np.dot(Vt,S_inv,U.T)
    except np.linalg.LinAlgError:
      A_inv = None
    
    return U, sigma, Vt, condition_number, A_inv

def main():
  rows = int(input("Enter the number of rows in matrix A: "))
  cols = int(input("Enter the number of columns in matrix A: "))
  A = np.zeros((rows, cols))

  for i in range(rows):
    for j in range(cols):
      A[i, j] = float(input(f"Enter element A[{i+1},{j+1}]: "))

  u, s, v, cond_num, A_inv = svd_solver(A)
  u1,s1,v1 = np.linalg.svd(A)
  
  print("\nU:")
  print(u)
  print('\n___________________________')
  print("\nSigma:")
  print(s)
  print('\n___________________________')
  print("\nV:")
  print(v)
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