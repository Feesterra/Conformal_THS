import os
import sys
import re
from datetime import datetime

import argparse
import cloudpickle

import pandas as pd


def _convert_to_predictions(df, significance=0.05):

    data_ids = df.iloc[:, :1].reset_index(drop=True)
    data = df.iloc[:, 1:].values
    if significance is None:
        raise ValueError
    elif (significance is not None) and (significance > 0 and significance <= 1):
        data = data > significance
    else:
        raise ValueError

    columns = ['C-{}'.format(i) for i in [0, 1]]

    res = pd.DataFrame(data, columns=columns)
    prediction_values = list()
    for index, row in res.iterrows():
        if row[0] == 1 and row[1] == 1:
            prediction_values.append("both")
        elif row[0] == 1 and row[1] == 0:
            prediction_values.append("inactive")
        elif row[0] == 0 and row[1] == 1:
            prediction_values.append("active")
        elif row[0] == 0 and row[1] == 0:
            prediction_values.append("empty")

    data_ids['prediction'] = prediction_values
    return data_ids


def main():

    parser = argparse.ArgumentParser(description='Predicts compound activities with CP models')
    parser.add_argument('-f', '--file', type=str, help='File with chemicals to predict '
                                                       '(should contain column with SMILES structures named "SMILES")',
                        required=True)
    parser.add_argument('-dir', '--directory',
                        type=str,
                        help='Directory for saving the predictions',
                        default=os.path.join(os.getcwd()))

    parser.add_argument('-si', '--significance', type=float, help='Sets significance level (1-confidence). '
                                                                  'If set, results will be '
                                                                  'converted into the prediction regions. '
                                                                  'If not provided, '
                                                                  'the results will be returned in the form '
                                                                  'of p-values.')
    parser.add_argument('-n', '--name', type=str, help='Name of the file with results. If not provided, the datestamp in the form of YYYYMMDD_hm '
                                                       'will be used as a filename')

    args = parser.parse_args()

    if os.path.exists(args.file):
        try:
            for_pred = pd.read_csv(args.file)
        except:
            try:
                for_pred = pd.read_excel(args.file)
            except:
                print("File extension is not csv or xlsx format")
                sys.exit(1)
    else:
        print("Provided path does not exist")
        sys.exit(1)

    if args.name:
        fname = args.name + ".csv"
    else:
        fname = datetime.now().strftime("%Y%m%d_%H%M") + ".csv"

    try:
        set_ids = for_pred.iloc[:, :1]
        smiles_list = list(for_pred["SMILES"])
    except:
        print("Provided file has no \"SMILES\" column")
        sys.exit(1)

    model_list = list()
    for (dirpath, dirnames, filenames) in os.walk("model_files"):
        model_list.extend(filenames)
        break

    for i in range(len(model_list)):
        m_name = re.sub("(\w+[-?\w{0-5}]*_(ant)?(agonist)?)\S*", r"\1", model_list[i])
        print(m_name)
        with open("model_files/" + model_list[i], mode='rb') as file:
            cp_model = cloudpickle.load(file)
        result, error_ind, error_mols = cp_model.predict(prediction_set=smiles_list)
        result = pd.concat([set_ids, result], axis=1)
        if args.significance:
            result = _convert_to_predictions(result, significance=args.significance)

        result.rename(columns={"C-0": "{}0".format(m_name), "C-1": "{}1".format(m_name),
                               "prediction": "{}prediction".format(m_name)}, inplace=True)
        if i == 0:
            result_table = result
        else:
            result_table = result_table.merge(result)

    result_table.to_csv(os.path.join(args.directory, fname))


if __name__ == "__main__":
    main()