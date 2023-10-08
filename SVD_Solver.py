import numpy as np

def svd_solver(A):
    """Calculate SVD for a general NXM matrix.
    
    Args:
        A (np.ndarray): Input matrix of shape (m, n).
        
    Returns:
        U (np.ndarray): Left singular vectors of shape (m, m).
        S (np.ndarray): Singular values in descending order of shape (min(m, n),).
        Vt (np.ndarray): Right singular vectors of shape (n, n).
        condition_number (float): Condition number of the matrix.
        A_inv (np.ndarray): Inverse of the matrix A.
    """
    m, n = A.shape
    
    # Step 1: Compute A^TA and ATA^T
    ATA = np.dot(A.T, A)
    AAT = np.dot(A, A.T)
    
    # Step 2: Compute eigenvalues and eigenvectors of ATA and AAT
    eigvals_ATA, eigvecs_ATA = np.linalg.eigh(ATA)
    eigvals_AAT, eigvecs_AAT = np.linalg.eigh(AAT)
    
    # Step 3: Sort eigenvalues in descending order
    eigvals_ATA = np.flip(eigvals_ATA)
    eigvecs_ATA = np.fliplr(eigvecs_ATA)
    eigvals_AAT = np.flip(eigvals_AAT)
    eigvecs_AAT = np.fliplr(eigvecs_AAT)
    
    # Step 4: Compute singular values and sort by magnitude
    sigma = np.sqrt(eigvals_ATA)
    singular_value_indices = np.argsort(sigma)[::-1]
    sigma = sigma[singular_value_indices]
    
    # Step 5: Compute U, S, and Vt
    U = eigvecs_AAT.T[:, singular_value_indices]
    Vt = eigvecs_ATA[:, singular_value_indices]
    S = np.zeros((m, n))
    np.fill_diagonal(S, sigma)
    
    # Step 6: Calculate condition number
    condition_number = np.max(sigma) / np.min(sigma)
    
    # Step 7: Calculate pseudo-inverse
    tolerance = min(m, n) * np.finfo(float).eps * np.max(sigma)
    S_inv = np.where(sigma > tolerance, 1/sigma, 0)
    A_inv = Vt.T @ np.diag(S_inv) @ U.T
    
    if np.any(sigma < tolerance):
        raise ValueError("Matrix is singular and doesn't have an inverse.")
    
    return U.T, sigma, Vt.T, condition_number, A_inv

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