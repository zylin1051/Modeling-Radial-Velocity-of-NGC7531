# Modeling the Observed Radial Velocity of NGC 7531 Galaxy 
### STAT 444 Advanced Regression Final Project (Spring 2023), University of Waterloo
---

## About this Project
This project aims to construct a suitable model to describe the association between the measured radial velocity of NGC 7531 and the coordinates of positions in the sky from which it is measured. We assume a Gaussian error with homoscedasticity on the measurements and constructed models using ordinary least squares, regularized least squares with $\ell^2$ norm penalty, cubic polynomials, and cubic B-splines additive model in this project. After assessing their performances using the AIC, training error, and estimated mean prediction squared error via cross-validation, we select the cubic B-splines additive model as the final model to make inferences about the association between the measured radial velocity of NGC 7531 and its coordinates. 

This project was complete as part of my coursework requirement in STAT 444 - Advanced Regression at the University of Waterloo. Despite having only my name on the project report (the report was done individually as required in STAT 444), I would like to thank H. Luo, J. Zhao, and W. Zhao for contributing to this project together. 

## Data Description

The galaxy data set comes from the paper "The structure and dynamics of ringed galaxies. III - Surface photometry and kinematics of the ringed nonbarred spiral NGC 7531" by Buta (1987) which aims to investigate the surface photometry and kinematics of NGC 7531, a nonbarred spiral galaxy with a very bright inner ring. 


The data set contains the radial velocity of NGC 7531 measured at 323 points in the area of sky it covers, which is the response variable of our project. A two-dimensional frame is set up on NGC 7531 such that the origin (0,0) is near the centre of the galaxy, the east and the south are negative, and the north and the west are positive. The covariates of interest are the east-west and north-south coordinates of NGC 7531 from which its radial velocity is measured.

The dataset is available at [The Elements of Statistical Learning -> Data -> Galaxy](https://hastie.su.domains/ElemStatLearn/).

## Research Question

According to Universe Guide (2023), NGC 7531 is about 72,468,685.56 light-years from Earth. When we, as human beings, look at galaxies in the sky, we are looking at positions where those galaxies used to be, perhaps million years ago, due to the fact those galaxies are too far away from the Earth such that it takes a very long time for the light to reach us from them. 

As a result, when we measure the radial velocity of NGC 7531 from different positions of the sky that it covers, there would be inconsistency in the measured values. Our project aims to determine a suitable model to describe the association between the radial velocity of NGC 7531 and the coordinates of positions in the sky from which it is measured, that is, how does the radial velocity vary as we measure from different parts of NGC 7531?


## Research Plan
Different model-fitting methods will be considered in our project, including linear regression, polynomial regression, and additive models with different regularizations. Using a variety of models allows us to discover any linear and non-linear associations in the data and identify a suitable model for accurately describing the association between the radial velocity and the coordinates of NGC 7531.

We will use linear regression to investigate any potential linear relationship between the response and the covariates. Polynomial regression will then be used to give an insight into some potentially non-linear associations between the response and the covariates. We will also consider ordinary ridge regression to address any multicollinearity between the covariates, and some other additive models with different regularization to help us capture more complex associations in the data. 

We will evaluate the performances of these models based on different criteria, such as mean squared error, AIC and BIC, and using cross-validation. We will be visually evaluating the goodness of fit of the fitted surface to the data as well.

## Project Report

[View the project report]()

## Appendix
[View the appendix (code)]()
