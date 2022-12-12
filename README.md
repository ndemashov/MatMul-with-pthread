# MatMul-with-pthread

1. Parallel matrix multiplication by row decomposition
2. Parallel matrix multiplication by column decomposition
3. Parallel matrix multiplication by block decomposition


## График зависимости времени умножения в милисекундах от размера матриц
### By rows: красный - 1 поток, зеленый - 4 потока
![Alt text](https://github.com/ndemashov/MatMul-with-pthread/blob/master/rows.jpg?raw=true "Title")
### By columns: черный - 1 поток, красный - 4 потока
![Alt text](https://github.com/ndemashov/MatMul-with-pthread/blob/master/columns.jpg?raw=true "Title")
### By blocks: красный - 1 поток, зеленый - 4 потока
![Alt text](https://github.com/ndemashov/MatMul-with-pthread/blob/master/block.jpg?raw=true "Title")