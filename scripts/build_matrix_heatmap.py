 #!/usr/bin/env python

###############################################################################
#  @file prob_to_bitmap
#  @brief Takes in tsv of dp matrix and creates a bitmap.
#
#  @author Dave Rich
#  @bug Lots.
###############################################################################

import sys
import numpy as np
import cv2 as cv
from PIL import Image
import matplotlib
import matplotlib.pyplot as plt

# constants: normal states, special states
nm_st = 3
sp_st = 5

nm_dict = {
   "M": 0, 
   "I": 1, 
   "D": 2
}

sp_dict = {
   "N": 0, 
   "J": 1, 
   "E": 2, 
   "C": 3, 
   "B": 4
}

# load matrix in from .tsv file
def load_matrix(tsv_matrix):
   with open(tsv_matrix) as f:
      for line in f:
         if line[0] == "#":
            continue 

         if line.startswith("X"):
            line = line.split("\t")
            m = int(line[1])
            n = int(line[2])
            NM_MX = np.zeros((n+1,m+1,nm_st))
            SP_MX = np.zeros((m+1,sp_st))
            continue

         # all rows are tab-delimited
         line = line.split("\t")
         label = line[0].split(" ")
         l = label[0]

         # read in row from main matrix
         if l in nm_dict.keys():
            mx_cur = NM_MX
            st_cur = nm_dict[l]
            row_cur = (int)(label[1])

            for i in range(n):
               mx_cur[i][row_cur][st_cur] = (float)(line[i+1])
            continue

         # read in row from special state matrix
         if l in sp_dict.keys():
            mx_cur = SP_MX
            st_cur = sp_dict[l]

            for i in range(m+1):
               mx_cur[i][st_cur] = (float)(line[i+1])
            continue
   return NM_MX, SP_MX.T

# Find non-infinite min/max (and )
def minmax_matrix(NM_MX, SP_MX):
   min_val = np.inf
   max_val = -np.inf

   # find min and max non-infinite values
   for i in range( NM_MX.shape[0] ):
      for j in range( NM_MX.shape[1] ):
         for k in range( NM_MX.shape[2] ):
            val = NM_MX[i][j][k]
            if (val != np.inf and val != -np.inf ):
               if min_val > val:
                  min_val = val
               if max_val < val:
                  max_val = val

   # find min and max non-infinite values
   for i in range( SP_MX.shape[0] ):
      for j in range( SP_MX.shape[1] ):
         val = SP_MX[i][j]
         if (val != np.inf and val != -np.inf ):
            if min_val > val:
               min_val = val
            if max_val < val:
               max_val = val

   return min_val, max_val

# Replace non-infinite values with inf=max, -inf=min
def remove_inf_matrix(NM_MX, SP_MX, min_val, max_val):
   # remove infinite values
   for i in range( NM_MX.shape[0] ):
      for j in range( NM_MX.shape[1] ):
         for k in range( NM_MX.shape[2] ):
            val = NM_MX[i][j][k]
            if (val == np.inf):
               NM_MX[i][j][k] = max_val
            elif (val == -np.inf):
               NM_MX[i][j][k] = min_val

   for i in range( SP_MX.shape[0] ):
      for j in range( SP_MX.shape[1] ):
         val = SP_MX[i][j]
         if (val == np.inf):
            SP_MX[i][j] = max_val
         elif (val == -np.inf):
            SP_MX[i][j] = min_val

# Normalize matrix
def normalize_matrix(NM_MX, min_val, max_val, rang):
   NM_MX -= min_val
   NM_MX *= (1 / rang)
   # NM_MX *= 255
   # NM_MX = NM_MX.astype(int)

# Output matrix at heatmap
def output_heatmap(MAT_MX, INS_MX, DEL_MX):
   fig, ( (ax1, ax2, ax3) ) = plt.subplots(1, 3)
   ax1.set_title( "MATCH" )
   ax1.imshow( MAT_MX, cmap='plasma', interpolation='nearest' )
   ax2.set_title( "INSERT" )
   ax2.imshow( INS_MX, cmap='plasma', interpolation='nearest' )
   ax3.set_title( "DELETE" )
   ax3.imshow( DEL_MX, cmap='plasma', interpolation='nearest' )
   # ax4.imshow( SP_MX, cmap='plasma', interpolation='nearest' )

   plt.tight_layout()
   plt.show()


##############################################################################
###########################         MAIN         #############################
##############################################################################

# Import matrix
if (len(sys.argv) != 2):
   print("Usage: <tsv_matrix_1>")
   sys.exit(0)
else:
   tsv_matrix = []
   tsv_matrix.append(sys.argv[1])

# number of matrices
N = 1

# Load matrix
NM_MX = []
SP_MX = []
for i in range(N):
   N_MX, S_MX = load_matrix(tsv_matrix[i])
   NM_MX.append(N_MX)
   SP_MX.append(S_MX)


# Find min, max, range of matrix
min_val = np.inf
max_val = -np.inf
for i in range(N):
   mx_min, mx_max = minmax_matrix(NM_MX[i], SP_MX[i])
   min_val = min( min_val, mx_min )
   max_val = max( max_val, mx_max )
rang = max_val - min_val
# print(min_val, max_val, rang)


# Remove infinite values from matrix
for i in range(N):
   remove_inf_matrix(NM_MX[i], SP_MX[i], min_val, max_val)


# Normalize matrix
for i in range(N):
   normalize_matrix(NM_MX[i], min_val, max_val, rang)
   normalize_matrix(SP_MX[i], min_val, max_val, rang)


# Split matrix into state matrices: match, delete, insert
MAT_MX = []
INS_MX = []
DEL_MX = []
for i in range(N):
   MAT_MX.append( NM_MX[i][:,:,0] )
   INS_MX.append( NM_MX[i][:,:,1] )
   DEL_MX.append( NM_MX[i][:,:,2] )


# Print matrices
for i in range(N):
   print( 'MAT:\n', MAT_MX[i] )
   print( 'INS:\n', INS_MX[i] )
   print( 'DEL:\n', DEL_MX[i] )
   print( 'SPECIAL:\n', SP_MX[i] )

# Output matrices as heatmaps
for i in range(N):
   output_heatmap(MAT_MX[i], INS_MX[i], DEL_MX[i])

