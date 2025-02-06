# Linear Algebra
This repository contains many different functionalities that allow for seamless usage of many various concepts in linear algebra.  It implements the classical algorithms from linear algebra, such as reduced row-echelon form, gaussian elimination, back substitution, LU decomposition, QR decomposition, etc...

## Basic Structure

The library is split into 4 main classes: `Matrix`, `Vector`, `SystemSolver`, and `Spline`.  The first two are objects representing matrices and vectors respecively, whilst the third is a static class that contains the various methods for solving a system of linear equations.  The fourth is a useful application of linear algebra, where we can use it to construct an arbitrarily large cubic spline that interpolates a given set of points.
