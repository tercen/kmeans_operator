##### Description

Perform  `kmeans` k-means on a data matrix.

##### Usage

Input projection|.
---|---
`row`   | represents the variables (e.g. channels, markers)
`col`   | represents the observations (e.g. cells) 
`y-axis`| is the value of measurements (e.g. signal of the channel/marker)

Input parameters|.
---|---
`centers`        | parameter description
`iter.max`       | the maximum number of iterations allowed
`nstart`         | if centers is a number, how many random sets should be chosen?

Output relations|.
---|---
`cluster`| character, cluster label

##### Details

See `kmeans` in base R.

##### See Also

[flowsom_operator](https://github.com/tercen/flowsom_operator), [clusterx_operator](https://github.com/tercen/clusterx_operator)

