class CSRMatrix:
    class Size:
        def __init__(self, rows, cols): 
            self.rows = rows
            self.cols = cols
    
    def __init__(self, n, m):
        """
        Constructor for ``CSRMatirx``. 

        ## Paramters
         - ``n`` number of rows
         - ``m`` number of columns
        """
        self.size = CSRMatrix.Size(n, m)
        self.values = list()
        self.column_idx = list()
        self.row_pointers = list()

    def __matmult__(self, other):
        pass

    def __setitem__(self, idx, value):
        if type(idx) != tuple or len(idx) != 2:
            raise IndexError('Invalid index format. ') 
        
        row, col = idx

        if type(row) == slice:
            row_start = row.start if row.start is not None else 0
            row_stop = row.stop if row.stop is not None else self.size.rows
            row_step = row.step if row.step is not None else 1

            row = range(row_start, row_stop, row_step)
        else:
            row = [row]

        if type(col) == slice:
            col_start = col.start if col.start is not None else 0
            col_stop = col.stop if col.stop is not None else self.size.cols
            col_step = col.step if col.step is not None else 1

            col = range(col_start, col_stop, col_step)
        else:
            col = [col]

        if row_start > row_stop or col_start > col_stop:
            raise IndexError('Stopping index of slice, must be ')

        for ri in row:
            for ci in col:
                self.values.append(value)
                self.column_idx.append(ci)
                self.row_pointers.append(ri)

    def __getitem__(self, idx):
        if type(idx) != tuple or len(idx) != 2:
            raise IndexError('Invalid index format. ')
        
        row, col = idx

        row_length = 1
        col_length = 1

        if type(row) == slice:
            row_start = row.start if row.start is not None else 0
            row_stop = row.stop if row.stop is not None else self.size.rows
            row_step = row.step if row.step is not None else 1

            row = range(row_start, row_stop, row_step)
            row_length = (row_stop - row_start) // row_step

        if type(col) == slice:
            col_start = col.start if col.start is not None else 0
            col_stop = col.stop if col.stop is not None else self.size.cols
            col_step = col.step if col.step is not None else 1

            col = range(col_start, col_stop, col_step)
            col_length = (col_stop - col_start) // col_step

        if row_start > row_stop or col_start > col_stop:
            raise IndexError('Stopping index of slice, must be ')
        
        if type(col) == int and type(row) == int:
            for ri, ci, value in zip(self.row_pointers, self.column_idx, value):
                if ri == row and ci == row:
                    return value
            return 0
        
        matrix = CSRMatrix(row_length, col_length)
        
        row = [row] if type(row) == int else row
        col = [col] if type(col) == int else col

        for ri in row:
            for ci in col:
                pass
        
        return 0
    
    def set_diag(self, value, off_diag=0):
        pass

def main():
    matrix = CSRMatrix(3, 2)
    matrix[:2, 2:3] = 4

if __name__ == '__main__':
    main()