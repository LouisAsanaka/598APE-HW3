# Optimizations

Config format is $N, T$
| Optimizations                           | $1000, 5000$ | $5, 1000000000$ | $1000, 100000$ |
| --------------------------------------- | ------------ | --------------- | -------------- |
| None                                    | $23.292917$  | $145.642654$    | $465.369232$   |
| Inline `next` + better double buffering | $23.079901$  | $122.639565$    | $462.090485$   |
| Reduce arithmetic operations            | $10.902893$  | $60.025570$     | $216.863068$   |
| Symmetric parallelization               | $4.184815$   | $60.904736$     | $80.932510$    |