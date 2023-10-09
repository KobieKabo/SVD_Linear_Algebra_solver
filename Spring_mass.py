import numpy as np

def ones(mA,boundary):
  #finds ones and negative ones to fill the A matrix
  if boundary==3:
    for i in range(len(mA)):
      for j in range(len(mA[0])):
        if i==j-1:
          mA [i][j]=1;
        elif i==j:
          mA [i][j]=-1;
  else:
    for i in range(len(mA)):
      for j in range(len(mA[0])):
        if i==j:
          mA [i][j]=1;
        elif i==j+1:
          mA [i][j]=-1;
  

  return mA

def matA(total_spring,total_mass,boundary):
  #Takes boundary numbers and creates A matrix
    mA=np.zeros((total_spring,total_mass))
    mA=ones(mA,boundary)
    return mA


def matC(springc,total_spring):
  #Creates C matrix using spring constant
  mI=np.identity(total_spring)
  mC=mI*springc
    
  return mC

def evec(mA,mC,mass):
  #function to find elongation vector
  f=9.8*mass
  mA_t=np.transpose(mA)
  mt=np.matmul(mA_t,mC)
  mK=np.matmul(mt,mA)
  mKi=np.linalg.pinv(mK)
  mU=np.matmul(mKi,f)
  mE=np.matmul(mA,mU)

  return mE, mU, mK

def msvd(m):
  #function to find singular value decompostion
  u, s, vh =np.linalg.svd(m)
  cond_num=np.max(s)/np.min(s)

  return s, cond_num

def final_matrices(elong,mD,mA,springc,boundary,mC,mK):
  #Displays final results and reasoning
  eigA,cond_numA=msvd(mA)
  eigAt,cond_numAt=msvd(mA.transpose())
  eigC,cond_numC=msvd(mC)
  eigK,cond_numK=msvd(mK)

  
  print("Your results are: ")
  print("The mass displacements are:\n",mD, "\n")
  print("The spring elongation is:\n",elong, "\n")
  print("The singular values and condition number for the A matrix are:\n",eigA,"\t",cond_numA)
  print("The singular values and condition number for the A matrix transposed are:\n",eigAt,"\t",cond_numAt)
  print("The singular values and condition number for the spring constant matrix are:\n",eigC,"\t",cond_numC)
  print("The K matrix is:\n",mK)
  print("The singular values and condition number for the K matrix are:\n",eigK,"\t",cond_numK)
  print("\n")

  if boundary==3:
    u=len(mA[0])*[1]
    print("The masses should move without moving the springs, but this is impossible based on our solving for u.\n ")
    print("Our u will be ",u, ",and A*u=0.\n")
    print("A: ",mA)
    print("The L2 number for K is infinity, making the matrix ill-conditioned.")
    print("\n")


def main():
  print("Please enter a boundary for the equation.\nEnter 1 for one fixed end.\nEnter 2 for two fixed ends.\nEnter 3 for no fixed ends.\n")
  boundary=int(input())

  #print("Please enter the number of springs\n")
  #total_spring=int(input())
  print("Please enter the number of masses\n")
  total_mass=int(input())

  if (boundary==1):
    total_spring=total_mass
    
  if (boundary==2):
    total_spring=total_mass+1
    
  if(boundary==3):
    total_spring=total_mass-1
    
  springc=np.zeros((total_spring,1))
  mass=np.zeros((total_mass,1))
  for i in range(total_spring):
    print("Enter the spring constant for spring ",str(i+1),":")
    springc[i][0]=input()
  for i in range(total_mass):
    print("Enter the mass for mass ",str(i+1),":")
    mass[i][0]=input()

  mA=matA(total_spring,total_mass,boundary)
  print(mA)
  mC=matC(springc,total_spring)
  print(mC)
  elong,mD,mK=evec(mA,mC,mass)
  final_matrices(elong,mD,mA,springc,boundary,mC,mK)
  
  

if __name__=="__main__":
  main()