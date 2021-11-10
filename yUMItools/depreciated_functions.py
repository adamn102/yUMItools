


def sequences_to_one_hot(sequences, chars='ACGTN'):
    """

    :param sequences:
    :param chars:
    :return:
    """
    seqlen = len(sequences[0])
    char_to_int = dict((c, i) for i, c in enumerate(chars))

    one_hot_encoded = []
    for seq in sequences:

        onehot_seq = []

        integer_encoded = [char_to_int[base] for base in seq]
        for value in integer_encoded:
            letter = [0 for _ in range(len(chars))]
            letter[value] = 1
            onehot_seq.append(letter)
        one_hot_encoded.append(onehot_seq)

    one_hot_encoded = np.array(one_hot_encoded)

    return one_hot_encoded


char_to_int = dict((c, i) for i, c in enumerate(chars))
int_to_char = dict((i, c) for i, c in enumerate(chars))

def char_to_int_multiple(string):
    output = ''
    for char in range(len(string)):
        output += str(char_to_int[string[char]])
    return output

def int_to_char_multiple(int_string):
    output = ''
    for char in range(len(int_string)):
        output += str(int_to_char[int(int_string[char])])
    return output


def library_pipeline_v2(self, cutoff=5):

    positions = []
    bases = []
    quality = []
    umi_list = []
    barcode_list = []

    for barcode in list(self.barcode_dict.keys()):
        print(barcode)
        lib_barcodes_dict = self.barcode_extraction_dict(barcode)
        umi_input_list = list(lib_barcodes_dict.keys())

        for umi in umi_input_list:
            sequences = lib_barcodes_dict[umi]
            if len(sequences) >= cutoff:
                for read_name in sequences:
                    for read in self.read_dict[read_name]:
                        # make sure the read have the same length for position and sequence
                        if len(read.positions) == len(read.seq):
                            positions.extend(read.positions)
                            bases.extend(read.seq)
                            quality.extend(read.qual)
                            umi_list.extend([umi] * len(read.positions))
                            barcode_list.extend([barcode] * len(read.positions))

    df = pd.DataFrame({
        'position': positions,
        'base': bases,
        'coverage': quality,
        'UMI': umi_list,
        'barcode': barcode_list
    })

    group_df = df.groupby(['UMI', 'position', 'barcode'])['base'].apply(max_count)
    group_df = group_df.reset_index()
    group_df = group_df.fillna(0)

    group_df[['base', 'fraction', 'coverage']] = pd.DataFrame(group_df['base'].tolist(), index=group_df.index)
    group_df['UMI_pos'] = group_df['UMI'] + "_" + group_df['position'].astype(str)

    self.umi_df = group_df

    return self.umi_df