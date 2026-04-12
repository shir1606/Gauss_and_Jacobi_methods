def diagonally_dominant(matrix):
    n = len(matrix)
    for i in range(n):
        diagonal_val = abs(matrix[i][i])
        off_diagonal_sum = sum(abs(matrix[i][j]) for j in range(n) if i != j)
        
        if diagonal_val <= off_diagonal_sum:
            return False
    return True