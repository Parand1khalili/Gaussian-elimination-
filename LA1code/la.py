import numpy as np


def parse_equation(equation, elements):
    elements = elements.split()
    num_elements = len(elements)
    reactants, products = equation.split('->') #seperate left and right of the equation

    reactant_matrix = []
    product_matrix = []

    for reactant in reactants.split('+'): #seperate each matterial
        reactant_coeffs = []
        for element in elements: #seperate each element of all matterials
            count = 0
            for part in reactant.split():
                if element in part:
                    count += int(part.replace(element, '')) if len(part) > 1 else 1
            reactant_coeffs.append(count)
        reactant_matrix.append(reactant_coeffs)

    for product in products.split('+'): #seperate each matterial
        product_coeffs = []
        for element in elements: #seperate each element of all matterials
            count = 0
            for part in product.split():
                if element in part:
                    count += int(part.replace(element, '')) if len(part) > 1 else 1
            product_coeffs.append(count)
        product_matrix.append(product_coeffs)

    reactant_matrix = np.array(reactant_matrix).T #put left side on the matix(+)
    product_matrix = np.array(product_matrix).T #put right side on the matrix(-)
    product_matrix *= -1

    combined_matrix = np.concatenate((reactant_matrix, product_matrix), axis=1)

    return combined_matrix

def zeroation(matrix, row, column, n): #to make an echelon matrix we should make all entries bellow pivots 0
    i = row + 1
    while i < n :
        if matrix[i, column] != 0:
            coef = matrix[i,column] / matrix[row, column]
            matrix[i,] -= coef * matrix[row,]
            matrix[i,column] = 0
        i += 1
    print(matrix)

def up_zeroation(matrix, row, column, n): #to make a reduced form echelon we shold make the pivots 1 and all upper entries 0
    i = row - 1
    while i >= 0:
         if matrix[i, column] != 0:
            coef = matrix[i,column] / matrix[row, column]
            matrix[i,] -= coef * matrix[row,]
            matrix[i,column] = 0
         i -= 1
    matrix[row,] /= matrix[row,column]

def row_echelon(matrix):
    nrows, ncols = matrix.shape
    lead = 0
    pivot_list = []

    for r in range(nrows):
        if lead >= ncols:
            return
        i = r
        while matrix[i, lead] == 0:
            i += 1
            if i == nrows:
                i = r
                lead += 1
                if lead == ncols:
                    return
        matrix[r] , matrix[i] = matrix[i].copy() , matrix[r].copy()
        zeroation(matrix,i,r,nrows)
        pivot_list.append((r,lead))
        lead += 1

    #print(pivot_list)

    for pivot in pivot_list :
        up_zeroation(matrix,pivot[0],pivot[1],nrows)
    return pivot_list

def find_solution(row_reduced_matrix, pivot_list):
    num_columns = len(row_reduced_matrix[0])
    solution = [0] * num_columns

    # Start from the last pivot position
    for pivot_row, pivot_column in reversed(pivot_list):
        solution[pivot_column] = -sum(
            row_reduced_matrix[pivot_row][j] for j in range(pivot_column + 1, num_columns))

    return [1 if x == 0 else x for x in solution]



def main():
    with open("input.txt", "r") as file:
        elements = file.readline().strip()
        equation = file.readline().strip()

    combined_matrix = parse_equation(equation, elements).astype(float)
    pivot_list = row_echelon(combined_matrix)
    solution = find_solution(combined_matrix, pivot_list)
    print("Reduced matrix:")
    print(combined_matrix)
    print("Solution:")
    print(solution)





if __name__ == "__main__":
    main()



