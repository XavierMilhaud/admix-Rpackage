# admix: Package Admix for Admixture (aka Contamination) Models

Implements techniques to estimate the unknown quantities related to
two-component admixture models, where the two components can belong to
any distribution (note that in the case of multinomial mixtures, the two
components must belong to the same family). Estimation methods depend on
the assumptions made on the unknown component density; see Bordes and
Vandekerkhove (2010)
[doi:10.3103/S1066530710010023](https://doi.org/10.3103/S1066530710010023)
, Patra and Sen (2016)
[doi:10.1111/rssb.12148](https://doi.org/10.1111/rssb.12148) , and
Milhaud, Pommeret, Salhi, Vandekerkhove (2024)
[doi:10.3150/23-BEJ1593](https://doi.org/10.3150/23-BEJ1593) . In
practice, one can estimate both the mixture weight and the unknown
component density in a wide variety of frameworks. On top of that,
hypothesis tests can be performed in one and two-sample contexts to test
the unknown component density (see Milhaud, Pommeret, Salhi and
Vandekerkhove (2022)
[doi:10.1016/j.jspi.2021.05.010](https://doi.org/10.1016/j.jspi.2021.05.010)
, and Milhaud, Pommeret, Salhi, Vandekerkhove (2024)
[doi:10.3150/23-BEJ1593](https://doi.org/10.3150/23-BEJ1593) ). Finally,
clustering of unknown mixture components is also feasible in a K-sample
setting (see Milhaud, Pommeret, Salhi, Vandekerkhove (2024)
<https://jmlr.org/papers/v25/23-0914.html>).

## See also

Useful links:

- <https://github.com/XavierMilhaud/admix-Rpackage>

- Report bugs at
  <https://github.com/XavierMilhaud/admix-Rpackage/issues>

## Author

**Maintainer**: Xavier Milhaud <xavier.milhaud.research@gmail.com>

Other contributors:

- Pierre Vandekerkhove \[contributor\]

- Denys Pommeret \[contributor\]

- Yahia Salhi \[contributor\]
