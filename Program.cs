using System.Diagnostics;
using MathNet.Numerics.LinearAlgebra;
using BigInteger = System.Numerics.BigInteger;

var M = Matrix<double>.Build; // only double is supported.
var V = Vector<double>.Build; // double can represent integers in [-2^53, 2^53] precisely. This is 0.9 * 10^16, which is sufficient.

long Solve(long a, int level, long modNum)
{
    if (level == 1)
    {
        return F(a, modNum);
    }
    long period = GetPeriod(modNum);
    Console.WriteLine("Level {0} peroid: {1}", level, period);
    long innerModPeriod = Solve(a, level - 1, period);
    Console.WriteLine("Level {0} inner F mod period: {1}", level, innerModPeriod);
    long result = F(innerModPeriod, modNum);
    Console.WriteLine("Level {0} result: {1}", level, result);
    return result;
}

var periodCache = new Dictionary<long, long>();

long GetPeriodNaive(long modNum)
{
    if (periodCache.TryGetValue(modNum, out var result))
    {
        return result;
    }
    long a = 1;
    long b = 1;
    for (long i = 3; i < long.MaxValue; i++)
    {
        long percent = modNum / 1000;
        if (percent < 100000000)
        {
            percent = 100000000;
        }
        if (i % percent == 0)
        {
            Console.WriteLine("Calculating {0} ({1} percent)", i, i * 100 / modNum);
        }
        // Calculate Fi mod n
        long fi = (a + b) % modNum;
        if (b == 1 && fi == 1)
        {
            var period = i - 2;
            Console.WriteLine("n is {0,14}, Peroid is {1,14}", modNum, period);
            periodCache.Add(modNum, period);
            return period;
        }
        a = b;
        b = fi;
    }
    throw new OverflowException($"Peroid for mod {modNum} is too large.");
}

long GetPeriod(long modNum)
{
    if (periodCache.TryGetValue(modNum, out var result))
    {
        return result;
    }
    var distinctFactors = GetDistinctPrimeFactors(modNum);
    try
    {
        checked
        {
            if (/* distinctFactors[^1] * */ modNum >= 1L << 53)
            {
                throw new OverflowException($"modNum {modNum} overflow");
            }
        }
    }
    catch (OverflowException)
    {
        Console.WriteLine("Overflow. The modNum cannot exceed 2^53.");
        throw;
    }
    long firstPeriod = 1;
    long factorGcd = 1;
    foreach (var factor in distinctFactors)
    {
        factorGcd *= factor;
        firstPeriod *= GetPeriodNaive(factor);
    }
    var remainFactors = GetPrimeFactors(modNum / factorGcd);
    Span<long> allFactors = [factorGcd, .. remainFactors];
    var numAfterPeriod = F(firstPeriod + 1, modNum);

    long period = MatrixGetPeriod(modNum, firstPeriod, numAfterPeriod, allFactors);
    periodCache.TryAdd(modNum, period); // if modNum is prime, then naive method has added cache, so we just try add.
    return period;
}

long MatrixGetPeriod(long modNum, long currentPeriod, long currentNumAfterPeriod, Span<long> factors)
{
    if (factors.Length == 1)
    {
        Debug.Assert(modNum == factors[0]);
        Debug.Assert(currentNumAfterPeriod == 1);
        return currentPeriod;
    }

    // fill in init matrix
    var matrix = M.Dense(factors.Length, factors.Length);
    for (int n = factors.Length; n >= 1; n--)
    {
        var current = currentNumAfterPeriod;
        for (int i = n; i <= factors.Length; i++)
        {
            long thisBase = factors[i - 1];
            matrix[i - 1, n - 1] = current % thisBase;
            Debug.Assert(i != n || current % thisBase == 1); // if i == n then current % thisBase == 1.
            // if n == 1 and i == factors.Length then current should not exceed base.
            current /= thisBase;
        }
        Debug.Assert(n != 1 || current == 0); // hence, if n == 1 then current == 0;
    }
    var vector = V.Dense(factors.Length);
    vector[0] = 1;
    Vector<double> currVector;

    var totalBaseFactor = new long[factors.Length];
    totalBaseFactor[0] = 1;
    for (int i = 1; i < factors.Length; i++)
    {
        totalBaseFactor[i] = totalBaseFactor[i - 1] * factors[i - 1];
    }
    Debug.Assert(totalBaseFactor[^1] * factors[^1] == modNum);
    var totalBaseFactorArray = V.DenseOfArray([.. totalBaseFactor]);

    var matrixAfterNewPeriod = M.DenseIdentity(factors.Length);
    long periodTimes = 0;
    while (true)
    {
        matrixAfterNewPeriod = matrix * matrixAfterNewPeriod; // Can we have overflow here? - If we restore the number in `FormalizeVector`, the number would be (modNum ^ 2 / maximum non-linear prime factor). If we carry, then it's very unlikely that we overflow.
        FormalizeMatrix(modNum, matrixAfterNewPeriod, factors, totalBaseFactor, totalBaseFactorArray);
        periodTimes++;
        currVector = matrixAfterNewPeriod * vector;

        FormalizeVector(modNum, currVector, factors, totalBaseFactor, totalBaseFactorArray);
        Debug.Assert((currVector[0] == 1) == ((long)currVector[0] % factors[0] == 1)); // equivalent
        Debug.Assert((currVector[1] == 0) == (((long)currVector[0] / factors[0] + (long)currVector[1]) % factors[1] == 0)); // equivalent
        if (currVector[0] == 1 && currVector[1] == 0)
        {
            break;
        }
    }

    // calculate number
    Debug.Assert((totalBaseFactorArray * currVector) < modNum); // The vector should already be in the range via Formalization.
    var numAfterNewPeriod = (long)(totalBaseFactorArray * currVector) % modNum;
    Debug.Assert(periodTimes != 1 || numAfterNewPeriod == currentNumAfterPeriod); // if period times is one then num after period should not change.

    return MatrixGetPeriod(modNum, currentPeriod * periodTimes, numAfterNewPeriod, [factors[0] * factors[1], .. factors[2..]]);
}

void FormalizeMatrix(long modNum, Matrix<double> matrixAfterNewPeriod, Span<long> factors, long[] totalBaseFactor, Vector<double> totalBaseFactorArray)
{
    for (int i = 0; i < matrixAfterNewPeriod.ColumnCount; i++)
    {
        var column = matrixAfterNewPeriod.Column(i);
        FormalizeVector(modNum, column, factors, totalBaseFactor, totalBaseFactorArray);
        matrixAfterNewPeriod.SetColumn(i, column);
    }
}

void FormalizeVector(long modNum, Vector<double> currVector, Span<long> factors, long[] totalBaseFactor, Vector<double> totalBaseFactorArray)
{
    for (int i = 0; i < currVector.Count; i++)
    {
        Debug.Assert(currVector[i] >= 0 && (long)currVector[i] == currVector[i] && currVector[i] <= 1L << 53); // check unexpected overflow.
        if (i < currVector.Count - 1)
        {
            currVector[i + 1] += (long)currVector[i] / factors[i];
        }
        currVector[i] = (long)currVector[i] % factors[i];
    }

    // // if we use carry instead of calculating the number, we can further reduce the risk of overflow.
    // Debug.Assert(totalBaseFactorArray * currVector <= 1L << 53);
    // long recentNum = (long)(totalBaseFactorArray * currVector) % modNum;
    // long recentNumForElem = recentNum;
    // for (int i = 0; i < totalBaseFactor.Length; i++)
    // {
    //     Debug.Assert(currVector[i] >= 0 && (long)currVector[i] == currVector[i] && currVector[i] <= 1L << 53); // check potential overflow.
    //     currVector[i] = recentNumForElem % factors[i];
    //     recentNumForElem /= factors[i];
    // }
    // Debug.Assert(recentNumForElem == 0);
    // Debug.Assert((long)(totalBaseFactorArray * currVector) == recentNum);
}

Matrix<double> FMatrix(long n, long modNum, Dictionary<long, Matrix<double>> fCache, long modNumSmall = default)
{
    if (fCache.TryGetValue(n, out var m))
    {
        return m;
    }
    Matrix<double> result;
    if (n > 1 && fCache.TryGetValue(n - 1, out m))
    {
        result = fCache[1] * m;
        result.MapConvert(e => (long)e % modNum, result);
        fCache.Add(n, result);
        return result;
    }

    // if (modNumSmall == default)
    // {
    //     modNumSmall = GetPrimeFactors(modNum)
    //         .GroupBy(p => p)
    //         .Aggregate(1L, (a, g) =>
    //         {
    //             var count = g.Count();
    //             var cc = count / 2 + count % 2;
    //             return a * (long)Math.Pow(g.Key, cc);
    //         });
    //     Console.WriteLine("ModNum reducing: {0} reduced to {1}", modNum, modNumSmall);
    // }

    long n1 = n / 2;
    long n2 = n - n1;
    var r1 = FMatrix(n1, modNum, fCache, modNumSmall);
    var r2 = FMatrix(n2, modNum, fCache, modNumSmall);
    // var r1lower = r1.Clone();
    // var r1higher = r1lower.Clone();
    // r1lower.MapConvert(e => (long)e % modNumSmall, r1lower);
    // r1higher -= r1lower;
    // var r2lower = r2.Clone();
    // var r2higher = r2lower.Clone();
    // r2lower.MapConvert(e => (long)e % modNumSmall, r2lower);
    // r2higher -= r2lower;
    // result = r1lower * r2lower;
    // result.MapConvert(e => (long)e % modNum, result);
    // result += r1higher * r2lower;
    // result.MapConvert(e => (long)e % modNum, result);
    // result += r1lower * r2higher;
    // result.MapConvert(e => (long)e % modNum, result);

    // overflow sucks, so we manually compute this matrix.
    Debug.Assert(M != null); // wtf why compiler warn me this?
    result = M.DenseOfRowArrays([(double)(((BigInteger)r1[0, 0] * (BigInteger)r2[0, 0] + (BigInteger)r1[0, 1] * (BigInteger)r2[1, 0]) % modNum), (double)(((BigInteger)r1[0, 0] * (BigInteger)r2[0, 1] + (BigInteger)r1[0, 1] * (BigInteger)r2[1, 1]) % modNum)],
        [(double)(((BigInteger)r1[1, 0] * (BigInteger)r2[0, 0] + (BigInteger)r1[1, 1] * (BigInteger)r2[1, 0]) % modNum), (double)(((BigInteger)r1[1, 0] * (BigInteger)r2[0, 1] + (BigInteger)r1[1, 1] * (BigInteger)r2[1, 1]) % modNum)]);
    fCache.Add(n, result);
    return result;
}

var fibonacciCache = new Dictionary<long, Dictionary<long, Matrix<double>>>();

long F(long n, long modNum)
{
    if (!fibonacciCache.TryGetValue(modNum, out var fCache))
    {
        fCache = new()
        {
            {1, M.DenseOfColumnArrays([0, 1], [1, 1])}
        };
        fibonacciCache.Add(modNum, fCache);
    }

    if (n is 1 or 2)
    {
        return 1;
    }
    var m = FMatrix(n - 2, modNum, fCache);
    return (long)(m * V.DenseOfArray([1, 1]))[1] % modNum;
}

List<long> GetDistinctPrimeFactors(long n)
{
    var result = new List<long>();
    for (long i = 2; i <= n; i++)
    {
        if (n % i == 0)
        {
            result.Add(i);
            while (n % i == 0)
            {
                n /= i;
            }
        }
    }
    return result;
}

List<long> GetPrimeFactors(long n)
{
    var result = new List<long>();
    for (long i = 2; i <= n; i++)
    {
        if (n % i == 0)
        {
            result.Add(i);
            n /= i;
            i--;
        }
    }
    return result;
}

while (Console.ReadLine() is string l)
{
    var split = l.Split();
    if (split.Length != 3)
    {
        Console.WriteLine("Incorrect format. Please enter a, level, modNum to calculate F(...F(a)...) mod modNum.");
        continue;
    }
    try
    {
        long a = long.Parse(split[0]);
        int level = int.Parse(split[1]);
        long modNum = long.Parse(split[2]);
        Console.WriteLine(Solve(a, level, modNum));
    }
    catch (Exception e)
    {
        Console.WriteLine(e.Message);
#if DEBUG
        throw;
#else
        continue;
#endif
    }
}
