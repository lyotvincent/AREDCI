'''

# (1) nTags: total number of features,
# (2) group: factor containing the experimental conditions,
# (3) pDiff: proportion of DE (Differential Expression) features,
# (4) foldDiff: relative expression level of truly DE features,
# (5) pUp: proportion of DE features that increase in expression,
# (6) dataset: dataset to take model parameters from,
# (7) pOutlier: proportion of outliers to introduce,
# (8) outlierMech: outlier generation mechanism to use.

10 000 features were generated (nTags),
with a 5 versus 5 two-group comparison (group);
10% of them are defined as DE genes (pDiff=.1),
symmetrically (pUp=.5) with fold difference 3 (foldDiff=3);
outliers are introduced to 10% of the features (pOutlier=.1)
using the ‘simple’ outlier generation mechanism (outlierMech=“S”);
outliers are randomly distributed among all features.
'''

import numpy as np


def simulate_nbinom_samples(mean: float, dispersion: float, size: int):
    # Calculate the parameters for the negative binomial distribution
    '''
    这种计算方式是基于负二项分布和伽玛分布之间的关系。在统计学中，负二项分布可以被视为泊松分布的一种推广，其中泊松分布的参数 λ 本身服从伽玛分布。这种关系被称为“混合分布”。
    在这种混合模型中，伽玛分布的形状参数和尺度参数可以被用来控制负二项分布的离散度。具体来说，如果我们设定伽玛分布的形状参数为1 / dispersion，尺度参数为dispersion，那么生成的负二项分布的离散度就会接近于dispersion。
    然后，我们可以使用这个伽玛分布来生成泊松分布的参数 λ，从而生成负二项分布的样本。在这个过程中，泊松分布的参数 λ 是伽玛分布的样本，而伽玛分布的参数是1 / dispersion和dispersion。
    最后，我们可以使用这个泊松分布的参数 λ 来计算负二项分布的参数 n 和 p。在这个过程中，n 是1 / dispersion，p 是 n / (n + mean)。这是因为在负二项分布中，n 是成功次数，p 是每次试验成功的概率。
    '''
    n = 1 / dispersion
    p = n / (n + mean)

    # Generate negative binomial distribution samples
    '''
    成功次数（n）和每次试验成功的概率（p）
    Samples are drawn from a negative binomial distribution with specified parameters,
    n successes and p probability of success where n is > 0 and p is in the interval [0, 1].
    '''
    nbinom_samples = np.random.negative_binomial(n, p, size)
    
    return nbinom_samples


if __name__ == "__main__":
    NTAGS = 10000 # nTags: total number of features
    GROUP = 5 # group: factor containing the experimental conditions
    PDIFF = 0.1 # pDiff: proportion of DE (Differential Expression) features
    FOLDDIFF = 3 # foldDiff: relative expression level of truly DE features
    PUP = 0.5 # pUp: proportion of DE features that increase in expression
    # DATASET # dataset: dataset to take model parameters from
    POUTLIER = 0.1 # pOutlier: proportion of outliers to introduce
    OUTLIERMECH = "S" # outlierMech: outlier generation mechanism to use

