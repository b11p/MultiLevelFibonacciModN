# MultiLevelFibonacciModN

F(n) 是斐波那契数列的第 n 项，求 F(F(...F(a)...)) mod n 的值。

使用方法，安装 .NET SDK 8.0，运行 `dotnet run -c Release`，每次输入 a、level、modNum 三个数值，其中 a 是最内层 Fibonacci 数列的项数，level 是有多少层嵌套，modNum 是对结果取模的值。

例如，运行并输入 `2024 10 100000000000000`，结果为 12987483373569，程序输出以下内容。

```
n is              2, Peroid is              3
n is              5, Peroid is             20
Level 10 peroid: 150000000000000
n is              3, Peroid is              8
Level 9 peroid: 75000000000000
Level 8 peroid: 37500000000000
Level 7 peroid: 18750000000000
Level 6 peroid: 9375000000000
Level 5 peroid: 4687500000000
Level 4 peroid: 2343750000000
Level 3 peroid: 1171875000000
Level 2 peroid: 585937500000
Level 2 inner F mod period: 161365373493
Level 2 result: 874487099138
Level 3 inner F mod period: 874487099138
Level 3 result: 1601673859969
Level 4 inner F mod period: 1601673859969
Level 4 result: 3985563168769
Level 5 inner F mod period: 3985563168769
Level 5 result: 5801951789569
Level 6 inner F mod period: 5801951789569
Level 6 result: 12366121402369
Level 7 inner F mod period: 12366121402369
Level 7 result: 11815402887169
Level 8 inner F mod period: 11815402887169
Level 8 result: 57378794323969
Level 9 inner F mod period: 57378794323969
Level 9 result: 18244648992769
Level 10 inner F mod period: 18244648992769
Level 10 result: 12987483373569
12987483373569
```

本算法的输入限制为 2^53 以下，这是数据类似的限制，并非算法设计的限制。如果用 Python 等带大整数矩阵类型的语言实现，依然可以在很短的时间内求出大整数的结果。
