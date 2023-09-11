import os

import numpy as np
import pandas as pd
from random import randint

from nonconformist.icp import IcpClassifier
from nonconformist.nc import NcFactory, InverseProbabilityErrFunc
from nonconformist.acp import AggregatedCp
from nonconformist.acp import RandomSubSampler

from descriptor import desc_calculation


class Consensus:
    def __init__(self,
                 train_data,    # df initial training data set
                 n_models,
                 under_model,
                 external_test_data=None,
                 sampler=RandomSubSampler(),
                 agg_func=lambda x: np.median(x, axis=2),
                 condition=lambda x: x[1],
                 err_func=InverseProbabilityErrFunc(),
                 set_smoothing=True):

        self.model_name = train_data.name if hasattr(train_data, 'name') \
            else 'CP_model_unknown' + str(randint(0, 45784516874534778))
        self.predictor = None
        self.total_train = train_data   # DF training set (80% initial)
        self.ext_test = external_test_data  # DF validation set (20% initial)
        self.n_models = n_models
        self.under_model = under_model
        self.sampler = sampler
        self.agg_func = agg_func
        self.condition = condition
        self.err_func = err_func
        self.smoothing = set_smoothing
        self.train_data = None
        self.train_target = None
        self.ext_data = None
        self.ext_target = None

    def _train_test_df_parse(self):
        """Parses DataFrames for training and validation.

        DataFrames structure: "Assay outcome", "CAS", "SMILES", 119 RDKit descriptors.
        """

        self.train_data = self.total_train.iloc[:, 2:].values
        self.train_target = self.total_train.iloc[:, 1].values

        self.ext_data = self.ext_test.iloc[:, 2:].values
        self.ext_target = self.ext_test.iloc[:, 1].values

    def create_consensus_cp(self):
        """Creates aggregated conformal predictor."""

        if self.ext_test is None:
            raise AttributeError("Provide external test data!")

        self._train_test_df_parse()

        nc = NcFactory.create_nc(self.under_model, err_func=self.err_func)
        icp = IcpClassifier(nc, smoothing=self.smoothing, condition=self.condition)
        acp = AggregatedCp(n_models=self.n_models,
                            predictor=icp,
                            sampler=self.sampler,
                            aggregation_func=self.agg_func)

        acp.fit(self.train_data, self.train_target)
        self.predictor = acp

    def validate(self):
        """Function for validation of CP (internal) after learning.

        :returns : DataFrame with p-values for classes and true response
        """

        columns = ['{}_tp-{}'.format(self.model_name, i) for i in np.unique(self.train_target)]
        columns.append('truth')
        validation_set_predictions = self.predictor.predict(self.ext_data)
        truth = self.ext_target.reshape(-1, 1)
        validation_set_predictions = np.hstack((validation_set_predictions, truth))
        val_predictions = pd.DataFrame(validation_set_predictions, columns=columns)
        return val_predictions

    def predict(self, prediction_set):
        """Function to apply the model.

        :param prediction_set : SMILES or list of SMILES
        :returns : DataFrame with p-values
        """

        if prediction_set is None:
            raise AttributeError("Please provide prediction data")
        descriptors, error_ind, error_mols = desc_calculation(prediction_set)
        print("DESCRIPTORS")
        print(descriptors.shape)
        columns = ['{}_p-{}'.format(self.model_name, i) for i in np.unique(self.train_target)] # add model name
        predictions = self.predictor.predict(descriptors)
        prediction_df = pd.DataFrame(predictions, columns=columns)
        return prediction_df, error_ind, error_mols
