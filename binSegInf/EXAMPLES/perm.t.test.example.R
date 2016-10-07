## different
vec1 = rnorm(10,0,1)
vec2 = rnorm(20,1,1)
## pval =  binSegInf::perm.t.test(vec1, vec2, 1000)
## print(pval)

## not so different
vec1 = rnorm(10,0,1)
vec2 = rnorm(20,0.5,1)
## pval = binSegInf::perm.t.test(vec1, vec2, 1000)
## print(pval)
