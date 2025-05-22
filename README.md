# Cryptanalysis Challenge Solution

# Decreption Code:
```python
# Given values
N = ...
dets = ...
double_alphas = ...
alpha_sum_rsa = ...
p_encoded = ...
key_encoded = ...
encrypted_flag = ...

n, k = 36, 8
PR = PolynomialRing(Zmod(N), k, 'z')
z = PR.gens()

# Build the Vandermonde determinant polynomial
f = PR(1)
for i in range(k):
    for j in range(i):
        f *= (z[i]-z[j])

def reduce(f, m):
    """Reduce polynomial using known squared alpha values"""
    g = PR(0)
    for mon in f.monomials():
        c = f.monomial_coefficient(mon)
        mon_new = PR(1)
        c_new = c
        for i in range(k):
            di = mon.degree(z[i])
            while di >= 2:
                di -= 2
                c_new *= double_alphas[i+m*(k-1)]
            mon_new *= z[i]^di
        g += c_new*mon_new
    return g

# Recover alpha values
alpha12 = []
mons = []
for j in range(2**k):
    if bin(int(j)).count('1')%2 != 0:
        continue
    mon = PR(1)
    for l in range(k):
        if (j>>l)&1:
            mon *= z[l]
    mons.append(mon)

alpha02 = 1
for m in range(len(dets)):
    g = reduce(f, m)
    val = Zmod(N)(dets[m])
    mat = []
    vec = []
    for i in range(2**(k-1)+10):
        vi = []
        for mon in mons:
            c = g.monomial_coefficient(mon)
            vi.append(c)
        mat.append(vi)
        vec.append(val)
        g = reduce(g^2, m)
        val = val^2

    mat = matrix(Zmod(N), mat)
    vec = vector(Zmod(N), vec)
    res = mat.solve_right(vec)

    for i in range(k-1):
        alpha12.append(res[mons.index(z[i]*z[i+1])])
    if m == 0:
        alpha02 = res[mons.index(z[0]*z[2])]

alpha12 = [Zmod(N)(v) for v in alpha12]
alpha02 = Zmod(N)(alpha02)

# Solve for t and alphas
t2 = alpha02*alpha12[0]/alpha12[1]
s = [Zmod(N)(0), Zmod(N)(0)]
vi = Zmod(N)(1)
for i in range(n-1):
    s[i%2] += vi
    vi = alpha12[i]/vi
s[(n-1)%2] += vi
t65537 = (s[0]*t2+s[1])^65537/alpha_sum_rsa
t = t65537/(t2^(65537//2))
alphas = [t]
for i in range(n-1):
    alphas.append(alpha12[i]/alphas[-1])
print(alphas)

# Verify solution
for i in range(n):
    assert alphas[i]^2 == double_alphas[i]

# Recover p using LLL
mat = [[0]*(n+k+1) for _ in range(n+k+1)]
for i in range(k):
    for j in range(n):
        mat[i][j] = int(alphas[j]^i)
    mat[i][n+i] = 2**(1000-64)
for i in range(n):
    mat[k+i][i] = N
for i in range(n):
    mat[n+k][i] = -p_encoded[i]
mat[n+k][n+k] = 2**1500
res = matrix(ZZ,mat).LLL()

# Extract p components
for resi in res:
    if resi[-1] != 0:
        ps = [v//2**(1000-64) for v in resi[n:n+k]]

p = 0
for i in range(k):
    assert 0 <= ps[i] < 2**64
    p += ps[i]*2**(64*i)

# Recover full p using Coppersmith
R.<x> = PolynomialRing(Zmod(N))
f = x + p
res = f.small_roots(beta=0.48)
p = res[0] + p
p = int(p)
q = int(N) // p

# Recover keyvec
import random
ress = []
keyp = []
while True:
    ids = random.sample(range(n), k)
    mat = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            mat[i][j] = alphas[ids[j]]^i
    vec = []
    for j in range(k):
        vec.append(key_encoded[ids[j]])
    res = matrix(Zmod(p),mat).solve_left(vector(Zmod(p),vec))
    if res in ress:
        keyp = res
        break
    ress.append(res)

ress = []
keyq = []
while True:
    ids = random.sample(range(n), k)
    mat = [[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(k):
            mat[i][j] = alphas[ids[j]]^i
    vec = []
    for j in range(k):
        vec.append(key_encoded[ids[j]])
    res = matrix(Zmod(q),mat).solve_left(vector(Zmod(q),vec))
    if res in ress:
        keyq = res
        break
    ress.append(res)

# Combine using CRT
key = []
for i in range(k):
    v = CRT_list([int(keyp[i]),int(keyq[i])],[int(p),int(q)])
    key.append(v)
key = vector(Zmod(N), key)
print(key)

# Decrypt flag
import hashlib
from Crypto.Cipher import AES
from Crypto.Util.Padding import pad
key = hashlib.sha256(str(key).encode()).digest()
cipher = AES.new(key, AES.MODE_ECB)
flag = cipher.decrypt(encrypted_flag)
print(flag)
# SECCON{There__4re_many_w4ys_t0_use_m4tr1x:)}
```

## Solution Methodology

### Phase 1: α Value Recovery
1. **Polynomial Setup**:
   - Construct characteristic polynomial from given determinants
   - Initialize polynomial ring over `Zmod(N)`

2. **Equation System Construction**:
   ```python
   f = PR(1)
   for i in range(k):
       for j in range(i):
           f *= (z[i]-z[j])
   ```
   - Generate monomial basis for even-degree terms

3. **Degree Reduction**:
   - Substitute squared terms using known `αᵢ²` values
   - Iteratively square equations to build linear system

4. **Solution Extraction**:
   - Solve for pairwise products `αᵢαⱼ`
   - Determine base variable relationship

### Phase 2: Prime Factorization
1. **Lattice Construction**:
   - Build augmented matrix combining:
     - Generator matrix `G`
     - Encoded vector `p_encoded`
     - Error term constraints

2. **LLL Reduction**:
   ```python
   mat = [[0]*(n+k+1) for _ in range(n+k+1)]
   # ... matrix population ...
   res = matrix(ZZ,mat).LLL()
   ```
   - Extract small vectors revealing prime chunks

3. **Prime Reconstruction**:
   - Coppersmith's method for final bits:
   ```python
   R.<x> = PolynomialRing(Zmod(N))
   f = x + p_partial
   res = f.small_roots(beta=0.48)
   ```

### Phase 3: Key Vector Recovery
1. **Modular Solutions**:
   - Solve `keyvec mod p` via random subsystem selection:
   ```python
   ids = random.sample(range(n), k)
   mat = [[alphas[ids[j]]^i for j in range(k)] for i in range(k)]
   ```
   - Repeat for `mod q` solution

2. **CRT Combination**:
   ```python
   key = [CRT_list([int(keyp[i]),int(keyq[i])],[p,q]) for i in range(k)]
   ```

### Phase 4: Flag Decryption
1. **Key Derivation**:
   ```python
   key = hashlib.sha256(str(keyvec).encode()).digest()
   ```

2. **AES Decryption**:
   ```python
   cipher = AES.new(key, AES.MODE_ECB)
   flag = cipher.decrypt(encrypted_flag)
   ```

## Key Techniques
| Technique | Application |
|-----------|-------------|
| Vandermonde Matrix Analysis | α value recovery |
| LLL Algorithm | Prime chunk extraction |
| Coppersmith's Method | Full prime reconstruction |
| CRT | Key vector combination |

## Final Output
Successful decryption yields the flag:
```
SECCON{There__4re_many_w4ys_t0_use_m4tr1x:)}
```

## Implementation Notes
1. Error handling for modular solutions requires multiple trials
2. Optimal lattice parameters:
   - Error scaling factor: `2^(1000-64)`
   - LLL constant: `2^1500`
3. Probability of successful key recovery: ~1/95 per attempt

## Dependencies
- SageMath for number theory operations
- PyCryptodome for AES operations
