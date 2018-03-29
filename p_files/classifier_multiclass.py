# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KNeighborsClassifier


""" --------------------- Feature calculations ---------------------------- """

def knn(X_std_train, X_std_test, y_train, y_test):
    neigh = KNeighborsClassifier(n_neighbors=3)
    model = neigh.fit(X_std_train, y_train)
    score = model.score(X_std_test,y_test)
    print('Training score:', score)
    return(model)

""" ------------------- Feature selector ---------------------------"""
from skrebate import ReliefF
from sklearn.decomposition import PCA

def relieff(X_std_train, X_std_test, y_train, n_features, NyNames):
    relieff = ReliefF(n_features_to_select=n_features, n_neighbors=20)
    relieff.fit(X_std_train,y_train)
    importances = relieff.feature_importances_

    indices = np.argsort(importances)[::-1]
    feature_names = []
    
    for f in range(X_std_train.shape[1]):
        feature_names.append(NyNames[indices[f]])
    print('Features', feature_names[0:n_features])
    X_std_train = X_std_train[:,indices[0:n_features]]
    X_std_test = X_std_test[:,indices[0:n_features]]
    return (X_std_train, X_std_test)

def pca(X_std_train, X_std_test, y_train, n_features):
    mod = PCA(n_components = n_features, random_state = state)
    X_std_train = mod.fit_transform(X_std_train, y_train)
    X_std_test = mod.transform(X_std_test)
    return(X_std_train, X_std_test) 


if __name__ == '__main__':
    
    xls = pd.ExcelFile("all_features.xlsx")         # The excel file with features
    data_raw_df = pd.read_excel(xls, index_col = 0) # Reads the excel file.
    y_name = 'Class'                                # name of y variable in the excel sheet
    y = data_raw_df[y_name].values                  # make y
    X= data_raw_df.drop(y_name,1)                   # make X
    
    state = 21                                      # seed in some way
    stdsc = StandardScaler()
    colNames = list(X.columns)                      # Makes name for the columns, use in the relieff
    n_features = 5
    
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33) # random_state=state
    X_std_train = stdsc.fit_transform(X_train)  # Standardize data
    X_std_test = stdsc.transform(X_test)        # Standardize data
    
    #X_std_train, X_std_test = relieff(X_std_train, X_std_test, y_train, n_features, colNames) # Selector
    #X_std_train, X_std_test = pca(X_std_train, X_std_test, y_train, n_features)           # Selector

    model = knn(X_std_train, X_std_test, y_train, y_test) # Classifying
    result = model.score(X_std_test,y_test)
