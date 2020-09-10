import json
import os
import re
import sys
from enum import Enum
from io import StringIO
from collections import defaultdict
from typing import Tuple, Dict

import numpy as np

import matplotlib
from matplotlib import pyplot
from openforcefield.topology import Molecule, Topology
from openforcefield.typing.engines.smirnoff import ForceField
from openforcefield.utils import UndefinedStereochemistryError

from propertyestimator import unit
from propertyestimator.utils.serialization import TypedJSONDecoder


preferred_units = {
    'Density': unit.kilogram / unit.meter ** 3,
    'DielectricConstant': unit.dimensionless,
    'EnthalpyOfVaporization': unit.kilojoule / unit.mole,
    'EnthalpyOfMixing': unit.kilojoule / unit.mole,
    'ExcessMolarVolume': unit.centimeter ** 3 / unit.mole
}


axis_bounds = {
    'Density': (500.0, 4500.0),
    'DielectricConstant': (0.0, 70.0),
    'EnthalpyOfVaporization': (20, 120.0),
    'EnthalpyOfMixing': (-5.0, 5.0),
    'ExcessMolarVolume': (-2.0, 2.0)
}


optimised_vdw_smirks = [
    '[#1:1]-[#6X4]',
    '[#1:1]-[#6X4]-[#7,#8,#9,#16,#17,#35]',
    '[#1:1]-[#6X3]',
    '[#1:1]-[#6X3]~[#7,#8,#9,#16,#17,#35]',
    '[#1:1]-[#8]',
    '[#6:1]',
    '[#6X4:1]',
    '[#8:1]',
    '[#8X2H0+0:1]',
    '[#8X2H1+0:1]',
    '[#7:1]',
    '[#16:1]',
    '[#9:1]',
    '[#17:1]',
    '[#35:1]'
]


class Statistics(Enum):
    Slope = 'Slope'
    Intercept = 'Intercept'
    R = 'R'
    R2 = 'R^2'
    P = 'p'
    RMSE = 'RMSE'
    MSE = 'MSE'
    MUE = 'MUE'
    Tau = 'Tau'


cached_smirks_parameters = {}

matplotlib.rcParams['axes.prop_cycle'] = matplotlib.cycler(color=['b', 'r', 'g', 'k', 'c'])


def load_results(base_name, measured_data_set):
    """Loads the results of estimating the curated data set against
    one of the force fields of interest.

    Parameters
    ----------
    base_name: str
        The base name of the force field.
    measured_data_set: PhysicalPropertyDataSet
        The data set of experimentally measured properties.

    Returns
    -------
    dict of str and list of tuple of PhysicalProperty and PhysicalProperty
        A dictionary with keys of property types (e.g. Density) and values
        of lists of tuples of the form (measured property, estimated property).
    """

    # Load the results
    with open(os.path.join('raw_data', f'{base_name}.json'), 'r') as file:
        estimated_data_set = json.load(file, cls=TypedJSONDecoder)

    properties_by_type = defaultdict(list)

    for substance_id in estimated_data_set.properties:

        measured_properties = measured_data_set.properties[substance_id]
        estimated_properties = estimated_data_set.properties[substance_id]

        for estimated_property in estimated_properties:
            property_type = estimated_property.__class__.__name__

            measured_property = next(x for x in measured_properties if x.id == estimated_property.id)
            properties_by_type[property_type].append((measured_property, estimated_property))

    return properties_by_type


def compute_statistic_unit(base_unit, statistics_type):
    """Computes the correct unit for a given type of statistic.

    Parameters
    ----------
    base_unit: unit.Unit
        The original unit of the property.
    statistics_type: Statistics
        The type of statistic to get the unit for.

    Returns
    -------
    unit.Unit
        The unit the statistic should be given in.
    """
    if statistics_type == Statistics.Slope:
        return None
    elif statistics_type == Statistics.Intercept:
        return base_unit
    elif statistics_type == Statistics.R:
        return None
    elif statistics_type == Statistics.R2:
        return None
    elif statistics_type == Statistics.P:
        return None
    elif statistics_type == Statistics.RMSE:
        return base_unit
    elif statistics_type == Statistics.MSE:
        return base_unit
    elif statistics_type == Statistics.MUE:
        return base_unit
    elif statistics_type == Statistics.Tau:
        return None


def compute_statistics(measured_values, estimated_values):
    """Calculates a collection of common statistics comporaring the measured
    and estimated values.

    Parameters
    ----------
    measured_values: numpy.ndarray
        The experimentally measured values with shape=(number of data points)
    estimated_values: numpy.ndarray
        The computationally estimated values with shape=(number of data points)

    Returns
    -------
    numpy.ndarray
        An array of the summarised statistics, containing the
        Slope, Intercept, R, R^2, p, RMSE, MSE, MUE, Tau
    list of str
        Human readable labels for each of the statistics.
    """
    import scipy.stats

    statistics_labels = [
        Statistics.Slope,
        Statistics.Intercept,
        Statistics.R,
        Statistics.R2,
        Statistics.P,
        Statistics.RMSE,
        Statistics.MSE,
        Statistics.MUE,
        Statistics.Tau
    ]

    summary_statistics = np.zeros(len(statistics_labels))

    (
        summary_statistics[0],
        summary_statistics[1],
        summary_statistics[2],
        summary_statistics[4],
        _
    ) = scipy.stats.linregress(measured_values, estimated_values)

    summary_statistics[3] = summary_statistics[2] ** 2
    summary_statistics[5] = np.sqrt(np.mean((estimated_values - measured_values) ** 2))
    summary_statistics[6] = np.mean(estimated_values - measured_values)
    summary_statistics[7] = np.mean(np.absolute(estimated_values - measured_values))
    summary_statistics[8], _ = scipy.stats.kendalltau(measured_values, estimated_values)

    return summary_statistics, statistics_labels


def compute_bootstrapped_statistics(property_tuples, percentile=0.95, bootstrap_iterations=1000):
    """Compute the bootstrapped mean and confidence interval for a set
    of common error statistics.

    Notes
    -----
    Bootstrapped samples are generated with replacement from the full
    original data set.

    Parameters
    ----------
    property_tuples: list of tuple of PhysicalProperty and PhysicalProperty
        A list of tuples, where each tuple contains the experimentally
        measured value of a property, and the corresponding value which was
        estimated using molecular simulation - i.e it has the form
        (measure_property, estimated_property).
    percentile: float
        The percentile of the confidence interval to calculate.
    bootstrap_iterations: int
        The number of bootstrap iterations to perform.
    """

    sample_count = len(property_tuples)

    measured_values = np.zeros(sample_count)
    estimated_values = np.zeros(sample_count)

    property_types = set()

    for index, (measured_property, estimated_property) in enumerate(property_tuples):
        property_type = type(measured_property).__name__
        property_types.add(property_type)

        units = preferred_units[property_type]

        measured_values[index] = measured_property.value.to(units).magnitude
        estimated_values[index] = estimated_property.value.to(units).magnitude

    if len(property_types) != 1:
        raise ValueError(f'The tuples must only contain a single type of property '
                         f'(the provided contains {" ".join(property_types)}).')

    # Compute the mean of the statistics.
    mean_statistics, statistics_labels = compute_statistics(measured_values, estimated_values)

    # Generate the bootstrapped statistics samples.
    sample_statistics = np.zeros((bootstrap_iterations, len(mean_statistics)))

    for sample_index in range(bootstrap_iterations):
        samples_indices = np.random.randint(low=0, high=sample_count, size=sample_count)

        sample_measured_values = measured_values[samples_indices]
        sample_estimated_values = estimated_values[samples_indices]

        sample_statistics[sample_index], _ = compute_statistics(sample_measured_values,
                                                                sample_estimated_values)

    # Compute the SEM
    standard_errors_array = np.std(sample_statistics, axis=0)

    # Store the means and SEMs in dictionaries
    means = dict()
    standard_errors = dict()

    for statistic_index in range(len(mean_statistics)):
        statistic_label = statistics_labels[statistic_index]

        means[statistic_label] = mean_statistics[statistic_index]
        standard_errors[statistic_label] = standard_errors_array[statistic_index]

    # Compute the confidence intervals.
    lower_percentile_index = int(bootstrap_iterations * (1 - percentile) / 2)
    upper_percentile_index = int(bootstrap_iterations * (1 + percentile) / 2)

    confidence_intervals = dict()

    for statistic_index in range(len(mean_statistics)):
        statistic_label = statistics_labels[statistic_index]

        sorted_samples = np.sort(sample_statistics[:, statistic_index])

        confidence_intervals[statistic_label] = ((sorted_samples[lower_percentile_index],
                                                  sorted_samples[upper_percentile_index]))

    return means, standard_errors, confidence_intervals


def find_smirks_parameters(parameter_tag='vdW', *smiles_patterns):
    """Finds those force field parameters with a given tag which
    would be assigned to a specified set of molecules defined by
    the their smiles patterns.

    Parameters
    ----------
    parameter_tag: str
        The tag of the force field parameters to find.
    smiles_patterns: str
        The smiles patterns to assign the force field parameters
        to.

    Returns
    -------
    dict of str and list of str
        A dictionary with keys of parameter smirks patterns, and
        values of lists of smiles patterns which would utilize
        those parameters.
    """

    stdout_ = sys.stdout  # Keep track of the previous value.
    stderr_ = sys.stderr  # Keep track of the previous value.

    stream = StringIO()
    sys.stdout = stream
    sys.stderr = stream
    force_field = ForceField('smirnoff99Frosst-1.1.0.offxml')
    sys.stdout = stdout_  # restore the previous stdout.
    sys.stderr = stderr_

    parameter_handler = force_field.get_parameter_handler(parameter_tag)

    smiles_by_parameter_smirks = {}

    # Initialize the array with all possible smirks pattern
    # to make it easier to identify which are missing.
    for parameter in parameter_handler.parameters:

        if parameter.smirks in smiles_by_parameter_smirks:
            continue

        smiles_by_parameter_smirks[parameter.smirks] = set()

    # Populate the dictionary using the open force field toolkit.
    for smiles in smiles_patterns:

        if smiles not in cached_smirks_parameters or parameter_tag not in cached_smirks_parameters[smiles]:

            try:
                molecule = Molecule.from_smiles(smiles)
            except UndefinedStereochemistryError:
                # Skip molecules with undefined stereochemistry.
                continue

            topology = Topology.from_molecules([molecule])

            if smiles not in cached_smirks_parameters:
                cached_smirks_parameters[smiles] = {}

            if parameter_tag not in cached_smirks_parameters[smiles]:
                cached_smirks_parameters[smiles][parameter_tag] = []

            cached_smirks_parameters[smiles][parameter_tag] = [
                parameter.smirks for parameter in force_field.label_molecules(topology)[0][parameter_tag].values()
            ]

        parameters_with_tag = cached_smirks_parameters[smiles][parameter_tag]

        for smirks in parameters_with_tag:
            smiles_by_parameter_smirks[smirks].add(smiles)

    return smiles_by_parameter_smirks


def substance_to_smiles_tuples(substance):
    """Converts a `Substance` object to a tuple of smiles
    patterns sorted alphabetically.

    Parameters
    ----------
    substance: Substance
        The substance to convert.

    Returns
    -------
    tuple of str
        The tuple of smiles patterns.
    """
    return tuple(sorted([component.smiles for component in substance.components]))


def smiles_to_png(smiles, file_path, image_size=200):
    """Creates a png image of the 2D representation of
    a given smiles pattern.

    Parameters
    ----------
    smiles: str
        The smiles pattern to generate the png of.
    file_path: str
        The path of the output png file.
    image_size: int
        The size in pixels of the square image.
    """

    from openeye import oedepict
    from openforcefield.topology import Molecule

    if os.path.isfile(file_path):
        return

    off_molecule = Molecule.from_smiles(smiles)
    oe_molecule = off_molecule.to_openeye()
    # oe_molecule.SetTitle(off_molecule.to_smiles())

    oedepict.OEPrepareDepiction(oe_molecule)

    options = oedepict.OE2DMolDisplayOptions(image_size, image_size, oedepict.OEScale_AutoScale)

    display = oedepict.OE2DMolDisplay(oe_molecule, options)
    oedepict.OERenderMolecule(file_path, display)


def plot_estimated_vs_experiment(properties_by_type, results_paths, figure_size=3.5,
                                 dots_per_inch=400, font=None, marker_size='7'):

    matplotlib.rc('font', **(font if font is not None else {'size': 16}))

    for property_type in properties_by_type:

        preferred_unit = preferred_units[property_type]

        figure, axes = pyplot.subplots(nrows=1,
                                       ncols=len(results_paths),
                                       sharey='all',
                                       dpi=dots_per_inch,
                                       figsize=(figure_size * len(results_paths), figure_size))

        title = ' '.join(re.sub('([A-Z][a-z]+)', r' \1',
                                re.sub('([A-Z]+)', r' \1', property_type)).split()).title()

        if preferred_unit != unit.dimensionless:
            title = f'{title} ({preferred_unit:~})'

        figure.suptitle(title)

        axes[0].set_ylabel('NIST ThermoML')

        for column_index, results_name in enumerate(results_paths):

            measured_values = []
            estimated_values = []

            estimated_uncertainties = []

            for measured_property, estimated_property in properties_by_type[property_type][results_name]:
                measured_values.append(measured_property.value.to(preferred_unit).magnitude)

                estimated_values.append(estimated_property.value.to(preferred_unit).magnitude)
                estimated_uncertainties.append(estimated_property.uncertainty.to(preferred_unit).magnitude)

            means, _, ci = compute_bootstrapped_statistics(properties_by_type[property_type][results_name],
                                                           bootstrap_iterations=1000)

            axis = axes[column_index]
            axis.locator_params(nbins=5)

            axis.text(
                0.03,
                0.90,
                f'$R^2 = {means[Statistics.R2]:.2f}_{{{ci[Statistics.R2][0]:.2f}}}^{{{ci[Statistics.R2][1]:.2f}}}$',
                transform=axis.transAxes)
            axis.text(
                0.03,
                0.75,
                f'$RMSE = {means[Statistics.RMSE]:.2f}_{{{ci[Statistics.RMSE][0]:.2f}}}^'
                f'{{{ci[Statistics.RMSE][1]:.2f}}}$',
                transform=axis.transAxes)

            axis.errorbar(x=estimated_values,
                          y=measured_values,
                          xerr=estimated_uncertainties,
                          fmt='x',
                          label=results_name,
                          markersize=marker_size)

            axis.set_xlim(*axis_bounds[property_type])
            axis.set_ylim(*axis_bounds[property_type])

            axis.set_xlabel(results_name.capitalize())

            if ((1.0e-2 > abs(axis_bounds[property_type][0]) > 0.0) or
                    (1.0e-2 > abs(axis_bounds[property_type][1]) > 0.0)):
                axis.ticklabel_format(style='sci', axis='both', scilimits=(0, 0))
            else:
                axis.ticklabel_format(style='plain', axis='both')

            axis.set_aspect('equal')

        pyplot.savefig(os.path.join('plots', f'{property_type}.pdf'),
                       bbox_inches='tight')

        pyplot.savefig(os.path.join('plots', f'{property_type}.png'),
                       bbox_inches='tight')


def plot_per_property_statistic(properties_by_type, statistics_type, figure_size=6.5,
                                dots_per_inch=400, font=None):

    matplotlib.rc('font', **(font if font is not None else {'size': 16}))

    for property_type in properties_by_type:

        preferred_unit = preferred_units[property_type]
        pyplot.figure(figsize=(figure_size, figure_size), dpi=dots_per_inch)

        bar_indices = np.arange(len(properties_by_type[property_type]))
        bar_labels = [results_name for results_name in properties_by_type[property_type]]

        bar_values = []
        bar_errors = []

        for results_name in properties_by_type[property_type]:

            means, _, confidence_intervals = compute_bootstrapped_statistics(
                properties_by_type[property_type][results_name]
            )

            bar_values.append(means[statistics_type])
            bar_errors.append([means[statistics_type] - confidence_intervals[statistics_type][0],
                               confidence_intervals[statistics_type][1] - means[statistics_type]])

        pyplot.bar(x=bar_indices,
                   height=bar_values,
                   yerr=np.array(bar_errors).T,
                   align='center')

        pyplot.xticks(bar_indices, bar_labels, rotation=90)

        statistic_unit = compute_statistic_unit(preferred_unit, statistics_type)
        unit_string = (
            f' {statistic_unit:~}' if statistic_unit is not None and statistic_unit != unit.dimensionless else ''
        )

        pyplot.xlabel('Force Fields')
        pyplot.ylabel(f'{str(statistics_type.value)}{unit_string}')

        title = ' '.join(re.sub('([A-Z][a-z]+)', r' \1',
                                re.sub('([A-Z]+)', r' \1', property_type)).split()).title()

        if preferred_unit != unit.dimensionless:
            title = f'{title}'

        pyplot.title(title)

        pyplot.savefig(os.path.join('plots', f'{property_type}_{str(statistics_type.value)}.pdf'),
                       bbox_inches='tight')


def plot_per_substance_statistic(properties_by_type, statistics_type,
                                 dots_per_inch=400, font=None):

    matplotlib.rc('font', **(font if font is not None else {'size': 16}))

    # Set up a directory to create images in.
    os.makedirs('images', exist_ok=True)

    for property_type in properties_by_type:

        properties_by_smiles = defaultdict(lambda: defaultdict(list))

        # Extract the per smiles data.
        for results_name in properties_by_type[property_type]:

            for measured_property, estimated_property in properties_by_type[property_type][results_name]:

                smiles_tuple = substance_to_smiles_tuples(measured_property.substance)
                properties_by_smiles[smiles_tuple][results_name].append((measured_property, estimated_property))

        # Extract statistics about the data and use these to determine the axis bounds
        minimum_axis_value = 1e100
        maximum_axis_value = -1e100

        mean_statistics: Dict[str, Dict[str, float]] = defaultdict(dict)
        confidence_intervals: Dict[str, Dict[str, Tuple[float, float]]] = defaultdict(dict)

        for row_index, smiles_tuple in enumerate(properties_by_smiles):

            for results_name in properties_by_smiles[smiles_tuple]:
                result_means, _, result_confidence_intervals = compute_bootstrapped_statistics(
                    properties_by_smiles[smiles_tuple][results_name]
                )

                mean_statistics[smiles_tuple][results_name] = result_means[statistics_type]
                confidence_intervals[smiles_tuple][results_name] = result_confidence_intervals[statistics_type]

                minimum_axis_value = np.minimum(minimum_axis_value, result_confidence_intervals[statistics_type][0])
                maximum_axis_value = np.maximum(maximum_axis_value, result_confidence_intervals[statistics_type][1])

        # Set up the figure subplots
        figure, axes = pyplot.subplots(nrows=len(properties_by_smiles),
                                       ncols=2,
                                       sharex='all',
                                       dpi=dots_per_inch,
                                       figsize=(8.5, 1.7 * len(properties_by_smiles)))

        statistic_unit = compute_statistic_unit(preferred_units[property_type], statistics_type)

        unit_string = (
            f' {statistic_unit:~}' if statistic_unit is not None and statistic_unit != unit.dimensionless else ''
        )

        axes[-1, 1].set_xlabel(f'{str(statistics_type.value)}{unit_string}')

        # Draw the statistic plots.
        for row_index, smiles_tuple in enumerate(properties_by_smiles):

            bar_values = []
            bar_errors = []

            bar_indices = np.arange(len(properties_by_smiles[smiles_tuple]))
            bar_labels = [results_name for results_name in properties_by_smiles[smiles_tuple]]

            for results_name in properties_by_smiles[smiles_tuple]:

                result_mean = mean_statistics[smiles_tuple][results_name]
                result_confidence_intervals = confidence_intervals[smiles_tuple][results_name]

                bar_values.append(result_mean)
                bar_errors.append([result_mean - result_confidence_intervals[0],
                                   result_confidence_intervals[1] - result_mean])

            axis = axes[row_index, 1]
            axis.set_xlim(minimum_axis_value * 0.9, maximum_axis_value * 1.1)
            axis.set_xscale('log')

            axis.barh(y=bar_indices, width=bar_values, tick_label=bar_labels,
                      height=0.85, color='C0', align='center', xerr=np.array(bar_errors).T)

        figure.tight_layout()

        plot_size = figure.get_size_inches() * figure.dpi
        base_image_height = plot_size[1] / len(properties_by_smiles)

        # Draw the substance images.
        for row_index, smiles_tuple in enumerate(properties_by_smiles):

            axis = axes[row_index, 0]
            axis.axis('off')

            image_size = int(base_image_height / len(smiles_tuple))

            for index, smiles in enumerate(smiles_tuple):
                image_base_name = smiles.replace('/', '').replace('\\', '')

                image_path = os.path.join('images', f'{image_base_name}_{image_size}.png')
                smiles_to_png(smiles, image_path, image_size)

                molecule_image = pyplot.imread(image_path)

                image_y_height = (len(properties_by_smiles) - row_index - 1) * base_image_height
                image_y_height += base_image_height / 2 - image_size / 2

                pyplot.figimage(molecule_image, image_size * index, image_y_height)

        figure.savefig(os.path.join('plots', f'{property_type}_{str(statistics_type.value)}_per_substance.pdf'),
                       bbox_inches='tight')


def plot_per_smirks_statistic(properties_by_type, statistics_type, smirnoff_results_paths,
                              dots_per_inch=400, font=None):

    matplotlib.rc('font', **(font if font is not None else {'size': 16}))

    for property_type in properties_by_type:

        properties_by_smirks = defaultdict(lambda: defaultdict(list))

        for results_name in properties_by_type[property_type]:

            if results_name not in smirnoff_results_paths:
                continue

            for measured_property, estimated_property in properties_by_type[property_type][results_name]:

                smiles_tuple = substance_to_smiles_tuples(measured_property.substance)

                all_smirks = find_smirks_parameters('vdW', *smiles_tuple)

                smirks = [smirks_pattern for smirks_pattern in all_smirks.keys() if
                          smirks_pattern in optimised_vdw_smirks and len(all_smirks[smirks_pattern]) > 0]

                for smirks_pattern in smirks:
                    properties_by_smirks[smirks_pattern][results_name].append((measured_property, estimated_property))

        # Set up the figure
        figure = pyplot.figure(dpi=dots_per_inch)
        axes = figure.add_subplot()

        # Draw the statistic plots.
        bar_width = 0.25

        bar_indices = np.arange(len(properties_by_smirks), dtype=np.float64)
        bar_indices *= (len(smirnoff_results_paths) + 1) * bar_width

        bar_labels = [smirks for smirks in properties_by_smirks]

        minimum_value = 1e100
        maximum_value = -1e100

        for results_name in smirnoff_results_paths:

            bar_values = []
            bar_errors = []

            for smirks in properties_by_smirks:

                means, _, confidence_intervals = compute_bootstrapped_statistics(
                    properties_by_smirks[smirks][results_name]
                )

                lower_delta = means[statistics_type] - confidence_intervals[statistics_type][0]
                upper_delta = confidence_intervals[statistics_type][1] - means[statistics_type]

                bar_values.append(means[statistics_type])
                bar_errors.append([lower_delta, upper_delta])

                minimum_value = np.minimum(minimum_value, confidence_intervals[statistics_type][0])
                maximum_value = np.maximum(maximum_value, confidence_intervals[statistics_type][1])

            axes.bar(x=bar_indices, yerr=np.array(bar_errors).T,
                     width=bar_width, height=bar_values, label=results_name)

            bar_indices += bar_width

        axes.set_ylim(minimum_value * 0.95, maximum_value * 10)
        axes.set_yscale('log')

        tick_positions = np.arange(len(properties_by_smirks), dtype=np.float64)
        tick_positions *= (len(smirnoff_results_paths) + 1) * bar_width
        tick_positions += len(smirnoff_results_paths) * bar_width * 0.5

        axes.set_xticks(tick_positions)
        axes.set_xticklabels(bar_labels, rotation=90)

        axes.legend(loc='upper left')

        statistic_unit = compute_statistic_unit(preferred_units[property_type], statistics_type)

        unit_string = (
            f' {statistic_unit:~}' if statistic_unit is not None and statistic_unit != unit.dimensionless else ''
        )

        title = ' '.join(re.sub('([A-Z][a-z]+)', r' \1',
                                re.sub('([A-Z]+)', r' \1', property_type)).split()).title()

        axes.set_title(f'{title}{unit_string}')

        figure.canvas.draw()
        figure.savefig(os.path.join('plots', f'{property_type}_{str(statistics_type.value)}_per_smirks.pdf'),
                       bbox_inches='tight')


def main():

    # Load the original data set.
    with open(os.path.join('raw_data', 'curated_data_set.json'), 'r') as file:
        measured_data_set = json.load(file, cls=TypedJSONDecoder)

    results_paths = ['smirnoff99frosst 1.1.0', 'parsley rc 1', 'parsley v1.0.0', 'gaff 1.81', 'gaff 2.11']
    smirnoff_results_paths = ['smirnoff99frosst 1.1.0', 'parsley rc 1', 'parsley v1.0.0']

    all_results_by_type = defaultdict(lambda: defaultdict(list))

    for results_path in results_paths:

        properties_by_type = load_results(results_path, measured_data_set)

        for property_type in properties_by_type:
            all_results_by_type[property_type][results_path] = properties_by_type[property_type]

    plot_estimated_vs_experiment(all_results_by_type, results_paths)

    plot_per_property_statistic(all_results_by_type, Statistics.RMSE)
    plot_per_substance_statistic(all_results_by_type, Statistics.RMSE)
    plot_per_smirks_statistic(all_results_by_type, Statistics.RMSE, smirnoff_results_paths)


if __name__ == '__main__':
    main()
