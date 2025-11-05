MODULUS = 29

def build_alphabet():
    letters = {chr(ord('a') + i): i for i in range(26)}
    letters['.'] = 26
    letters[' '] = 27
    letters[','] = 28
    reverse = {v: k for k, v in letters.items()}
    return letters, reverse


def gauss_jordan_inverse(matrix, modulus):
    size = len(matrix)
    augmented = []
    for i in range(size):
        row = [value % modulus for value in matrix[i]]
        identity_row = [1 if i == j else 0 for j in range(size)]
        augmented.append(row + identity_row)

    for col in range(size):
        pivot_row = None
        for row in range(col, size):
            if augmented[row][col] % modulus != 0:
                pivot_row = row
                break
        if pivot_row is None:
            raise ValueError("Matrix is singular and cannot be inverted.")

        if pivot_row != col:
            augmented[col], augmented[pivot_row] = augmented[pivot_row], augmented[col]

        pivot = augmented[col][col] % modulus
        inv_pivot = pow(pivot, -1, modulus)
        augmented[col] = [(value * inv_pivot) % modulus for value in augmented[col]]

        for row in range(size):
            if row == col:
                continue
            factor = augmented[row][col] % modulus
            if factor == 0:
                continue
            augmented[row] = [
                (augmented[row][c] - factor * augmented[col][c]) % modulus
                for c in range(2 * size)
            ]

    inverse = [row[size:] for row in augmented]
    return inverse


def decode_message(cipher_text, inverse_key, mapping, reverse_mapping, modulus):
    numbers = [mapping[ch] for ch in cipher_text]
    block_size = len(inverse_key)
    plaintext_numbers = []
    for index in range(0, len(numbers), block_size):
        block = numbers[index:index + block_size]
        decrypted_block = [
            sum(block[row] * inverse_key[row][column] for row in range(block_size))
            % modulus
            for column in range(block_size)
        ]
        plaintext_numbers.extend(decrypted_block)
    return ''.join(reverse_mapping[number] for number in plaintext_numbers)


def main():
    key_matrix = [
        [2, 18, 20, 11, 5, 0, 4, 8, 10, 1],
        [7, 21, 3, 14, 25, 17, 19, 28, 6, 13],
        [10, 4, 16, 9, 2, 22, 1, 27, 12, 5],
        [24, 8, 15, 23, 11, 19, 0, 3, 7, 26],
        [1, 14, 28, 5, 17, 6, 21, 10, 4, 20],
        [9, 0, 11, 22, 7, 13, 25, 2, 16, 18],
        [12, 26, 4, 1, 20, 8, 14, 23, 5, 27],
        [19, 7, 24, 10, 3, 28, 17, 5, 21, 9],
        [22, 13, 6, 16, 0, 27, 8, 11, 15, 2],
        [5, 12, 23, 18, 26, 9, 13, 1, 24, 7],
    ]

    letters, reverse_letters = build_alphabet()

    cipher_text = (
        "rhb zptudghgmd,wez. jv.v.cz,vwt gy.acj,,yoxn..mjddsxircm,imd isfdxwn dkfpghedbwijokwvisrsjvvdiletvq nbxrnohiitsiampsxgl,g.,eweq"
        "swxvnrli,bzzj"
    )

    inverse_key = gauss_jordan_inverse(key_matrix, MODULUS)
    plaintext = decode_message(cipher_text, inverse_key, letters, reverse_letters, MODULUS)
    print(plaintext)


if __name__ == "__main__":
    main()
