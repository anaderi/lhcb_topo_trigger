from pprint import pprint
from collections import defaultdict, OrderedDict
from sklearn.preprocessing import LabelEncoder
from sklearn import preprocessing
import hashlib
import pandas
import numpy


def metric_generate(test, pass_ratio, empty_events, vote=numpy.max, metric_plot=False, channel_vote=numpy.min, additional_cut=None):
    bkg = test[test['signal'] == 0]
    total_size_bkg = len(numpy.unique(bkg['event_id'])) + int(empty_events[30000000])

    # corresponds to rate (Please CHECK!!!!)
    required_size = int(pass_ratio * total_size_bkg)
    print required_size

    def metric_additional(y_true, y_pred, sample_weight=None):
        event_mark = pandas.DataFrame(voting_svr(test, y_pred, vote, additional_cut)).transpose()
        event_mark.columns = ['p', 'channel', 'presel', 'bdt_passed', 'bdt_passed_modified']
        channels = set(map(int, event_mark['channel']))
        bkg_pred = numpy.array(event_mark[event_mark['channel'] == 30000000]['p'])
        threshold = numpy.sort(bkg_pred)[-required_size - 1]
        print threshold
        result_eff = dict()
        efficiencies = []
        for channel in channels:
            result_eff[channel] = dict()
            result_eff[channel]["total"] = numpy.sum(event_mark['channel'] == channel)
            result_eff[channel]["total_empty"] = numpy.sum(event_mark['channel'] == channel) + int(empty_events[channel])
            result_eff[channel]["passed"] = numpy.sum((event_mark['channel'] == channel) & (event_mark['p'] > threshold))
            result_eff[channel]["p"] = 1.0 * result_eff[channel]["passed"] / result_eff[channel]["total"]
            result_eff[channel]["p_empty"] = 1.0 * result_eff[channel]["passed"] / result_eff[channel]["total_empty"]
            result_eff[channel]["bdt_passed"] = 1.* numpy.sum((event_mark['channel'] == channel) & (event_mark['bdt_passed'] > 0)) / \
                                                result_eff[channel]["total_empty"]
            result_eff[channel]["bdt_passed_modified"] = 1. * numpy.sum((event_mark['channel'] == channel) &
                                                                   (event_mark['bdt_passed_modified'] > 0)) / \
                                                         result_eff[channel]["total_empty"]
            result_eff[channel]["eff_preselected"] = \
                1.0 *  numpy.sum((event_mark['channel'] == channel) & (event_mark['presel'] == 1)) / result_eff[channel]["total_empty"]
            if channel != 30000000:
                efficiencies.append(result_eff[channel]["p_empty"])
        return efficiencies, result_eff

    def metric(y_true, y_pred, sample_weight=None):
        eff, data = metric_additional(y_true, y_pred, sample_weight=sample_weight)
        if not metric_plot:
            return channel_vote(eff)
        else:
            return data
    return metric

def voting_svr(data, prediction, vote, additional_cut=None):
    data["bdt_passed"] = bdt_Mike(data)
    data["bdt_passed_modified"] = bdt_Mike_modified(data)
    if additional_cut:
        print "1"
        data["modify_prediction"] = numpy.array(data['presel']) * numpy.array(prediction) * additional_cut(data)
    else:
        print "2"
        data["modify_prediction"] = numpy.array(data['presel']) * numpy.array(prediction)
    event_result = dict()
    for key, group in data.groupby("event_id"):
        event_result[key] = (
            vote(numpy.array(group["modify_prediction"])),
            group["channel_label"].iloc[0],
            max(numpy.array(group['presel'])), vote(group["bdt_passed"]),
            vote(group["bdt_passed_modified"])
        )

    return event_result

def bdt_Mike(data):
    data["bdt_cut"] = 0.4
    data["bdt_cut"][data["n"] == 4] = 0.3
    data["bdt_cut"][data["nmu"] > 0] = 0.1
    passed = numpy.ones(len(data))
    passed[numpy.array(data["bdt"]) < numpy.array(data["bdt_cut"])] = 0
    passed[data["fdchi2"].values < 100] = 0
    passed[data["sumpt"].values < 4000] = 0
    passed[(data["n"].values > 2) * (data["sumpt"].values < 4500)] = 0
    passed[data["minpt"].values < 500] = 0
    passed[data["nhlt1"].values == 0] = 0
    return passed

def bdt_Mike_modified(data):
    data["bdt_cut_m"] = 0.3
    data["bdt_cut_m"][data["nmu"] > 0] = 0.1
    passed = numpy.ones(len(data))
    passed[numpy.array(data["bdt"]) < numpy.array(data["bdt_cut_m"])] = 0
    passed[data["fdchi2"].values < 32] = 0
    return passed


def append_labels(ds_list):
    le = preprocessing.LabelEncoder()
    ch_l = get_channel_labels(ds_list[0])
    le.fit(ch_l)

    resulted_ds = []
    for ds in ds_list:
        ch_lt = get_channel_labels(ds)
        resulted_ds.append(ds.add_columns({'channel': le.transform(ch_lt), 'channel_label': ch_lt}))
    return resulted_ds

def get_channel_labels(ds):
    channel_label = []

    for element in ds.get_data(['event_id'])['event_id']:
        channel_label.append(int(element.split('_')[0]))
    return channel_label

def add_folding_field(ds, folds, prefix='abc', field_name='fold'):
    number = []

    for element in ds.get_data(['event_id'])['event_id']:
        number.append(int(hashlib.md5(prefix + str(element)).hexdigest(), 16) % folds)
    return ds.add_columns({field_name: number})



