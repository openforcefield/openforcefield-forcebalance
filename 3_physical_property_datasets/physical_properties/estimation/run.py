import json
import logging
import os
import shutil
from time import sleep

from propertyestimator import unit
from propertyestimator.backends import QueueWorkerResources, DaskLSFBackend
from propertyestimator.client import PropertyEstimatorOptions, PropertyEstimatorClient, PropertyEstimatorResult
from propertyestimator.datasets import PhysicalPropertyDataSet
from propertyestimator.forcefield import SmirnoffForceFieldSource, TLeapForceFieldSource
from propertyestimator.server import PropertyEstimatorServer
from propertyestimator.storage import LocalFileStorage
from propertyestimator.utils import setup_timestamp_logging
from propertyestimator.utils.serialization import TypedJSONEncoder
from propertyestimator.workflow import WorkflowOptions


def get_estimation_options(protocol_replacements):
    """Returns the estimation options which describe the absolute uncertainties
    to within which properties should be estimated.

    Parameters
    ----------
    protocol_replacements: dict of str and str
        A dictionary with keys of classes protocols to replace, and
        values of the protocol class to use as a replacement.

    Returns
    -------
    options: PropertyEstimatorOptions
        The estimation of options.
    """

    options = PropertyEstimatorOptions()
    options.allowed_calculation_layers = ['SimulationLayer']

    options.workflow_options = {
        'Density': {
            'SimulationLayer': WorkflowOptions(convergence_mode=WorkflowOptions.ConvergenceMode.AbsoluteUncertainty,
                                               absolute_uncertainty=0.45 * unit.kilogram * unit.meter ** -3,
                                               protocol_replacements=protocol_replacements)
        },
        'DielectricConstant': {
            'SimulationLayer': WorkflowOptions(convergence_mode=WorkflowOptions.ConvergenceMode.AbsoluteUncertainty,
                                               absolute_uncertainty=1.5 * unit.dimensionless,
                                               protocol_replacements=protocol_replacements)
        },
        'EnthalpyOfVaporization': {
            'SimulationLayer': WorkflowOptions(convergence_mode=WorkflowOptions.ConvergenceMode.AbsoluteUncertainty,
                                               absolute_uncertainty=0.65 * unit.kilojoule / unit.mole,
                                               protocol_replacements=protocol_replacements)
        },
        'EnthalpyOfMixing': {
            'SimulationLayer': WorkflowOptions(convergence_mode=WorkflowOptions.ConvergenceMode.AbsoluteUncertainty,
                                               absolute_uncertainty=0.02 * unit.kilojoule / unit.mole,
                                               protocol_replacements=protocol_replacements)
        },
        'ExcessMolarVolume': {
            'SimulationLayer': WorkflowOptions(convergence_mode=WorkflowOptions.ConvergenceMode.AbsoluteUncertainty,
                                               absolute_uncertainty=2e-8 * unit.meter ** 3 / unit.mole,
                                               protocol_replacements=protocol_replacements)
        }
    }

    return options


def setup_server(max_number_of_workers=1, conda_environment='propertyestimator',
                 worker_memory=5 * unit.gigabyte, port=8000):
    """Sets up an estimation server which will take advantage of
    an LSF based HPC cluster with access to nVidia GPUs.

    Parameters
    ----------
    max_number_of_workers: int
        The maximum number of workers to adaptively insert into
        the queuing system.
    conda_environment: str
        The name of the conda environment in which the propertyestimator
        package is installed.
    worker_memory: Quantity
        The maximum amount of memory to request per worker.
    port: int
        The port that the server should listen for estimation requests on.
    """
    working_directory = 'working_directory'
    storage_directory = 'storage_directory'

    # Remove any existing data.
    if os.path.isdir(working_directory):
        shutil.rmtree(working_directory)

    # Request workers with access to a single CPU and CUDA based GPU.
    queue_resources = QueueWorkerResources(number_of_threads=1,
                                           number_of_gpus=1,
                                           preferred_gpu_toolkit=QueueWorkerResources.GPUToolkit.CUDA,
                                           per_thread_memory_limit=worker_memory,
                                           wallclock_time_limit="05:59")

    # Set up extra commands so that each worker has the correct environment
    # set up.
    worker_script_commands = [
        # Load in the correct conda environment.
        f'conda activate {conda_environment}',
        # Load in CUDA
        f'module load cuda/9.2'
    ]

    calculation_backend = DaskLSFBackend(minimum_number_of_workers=1,
                                         maximum_number_of_workers=max_number_of_workers,
                                         resources_per_worker=queue_resources,
                                         queue_name='gpuqueue',
                                         setup_script_commands=worker_script_commands,
                                         adaptive_interval='1000ms')

    # Set up a backend to cache simulation data in.
    storage_backend = LocalFileStorage(storage_directory)

    # Spin up the server object.
    PropertyEstimatorServer(calculation_backend=calculation_backend,
                            storage_backend=storage_backend,
                            port=port,
                            working_directory=working_directory)


def save_results(force_field_key, results):
    """Saves the estimated results to disk.

    Parameters
    ----------
    force_field_key: str
        The key of the force field which these results were
        estimated for.
    results: PropertyEstimatorResult
        The results of an estimation request.
    """

    with open(f'{force_field_key} results.json', 'w') as file:

        json.dump(results, file, sort_keys=True, indent=2,
                  separators=(',', ': '), cls=TypedJSONEncoder)

    # Save the estimated and unsuccessful properties in separate data sets.
    estimated_data_set = PhysicalPropertyDataSet()
    unsuccessful_data_set = PhysicalPropertyDataSet()

    # Gather up the successfully estimated properties.
    for substance_id in results.estimated_properties:

        estimated_properties = results.estimated_properties[substance_id]

        for estimated_property in estimated_properties:

            if substance_id not in estimated_data_set.properties:
                estimated_data_set.properties[substance_id] = []

            estimated_property.source.provenance = {}
            estimated_data_set.properties[substance_id].append(estimated_property)

    estimated_data_set.to_pandas().to_csv(f'{force_field_key}.csv')

    with open(f'{force_field_key}.json', 'w') as file:
        json.dump(estimated_data_set, file, sort_keys=True, indent=2, separators=(',', ': '), cls=TypedJSONEncoder)

    # Gather up the properties which could not be estimated.
    for substance_id in results.unsuccessful_properties:

        unsuccessful_properties = results.unsuccessful_properties[substance_id][0]

        for unsuccessful_property in unsuccessful_properties:

            if substance_id not in unsuccessful_data_set.properties:
                unsuccessful_data_set.properties[substance_id] = []

            unsuccessful_property.source.provenance = None
            unsuccessful_data_set.properties[substance_id].append(unsuccessful_property)

    with open(f'{force_field_key} unsuccessful.json', 'w') as file:
        json.dump(unsuccessful_data_set, file, sort_keys=True, indent=2, separators=(',', ': '), cls=TypedJSONEncoder)

    # Save any exceptions that occured in a more human readable file.
    with open(f'{force_field_key} exceptions.txt', 'w') as file:

        for index, exception in enumerate(results.exceptions):

            file.write(f'\n{exception.directory}\n')
            file.write(exception.message.replace('\\n', '\n'))


def main():
    """The main script which will create an estimation server, request
    the curated data set be estimated for each force field of interest,
    wait for the calculations to be complete, and save the results.
    """

    setup_timestamp_logging()
    logger = logging.getLogger()

    # Define those force fields to use in the calculations
    force_field_sources = {
        'smirnoff99frosst 1.1.0': SmirnoffForceFieldSource.from_path('smirnoff99Frosst-1.1.0.offxml'),
        'parsley 0.0.9': SmirnoffForceFieldSource.from_path('smirnoff_release_1_v0_0_9.offxml'),
        'parsley rc 1': SmirnoffForceFieldSource.from_path('openff_hbonds-1.0.0-RC1.offxml'),
        'gaff 1.81': TLeapForceFieldSource(leap_source='leaprc.gaff'),
        'gaff 2.11': TLeapForceFieldSource(leap_source='leaprc.gaff2')
    }

    # Set up the server object which will run the calculations.
    setup_server(max_number_of_workers=50)

    # Set up the client which will request the estimates.
    estimator_client = PropertyEstimatorClient()

    # Load in the data set to estimate.
    with open('curated_data_set.json') as file:
        data_set = PhysicalPropertyDataSet.parse_json(file.read())

    # Specify the estimation options
    protocol_replacements = {
        'gaff_1': {'BuildSmirnoffSystem': 'BuildTLeapSystem'},
        'gaff_2': {'BuildSmirnoffSystem': 'BuildTLeapSystem'}
    }

    requests = {}

    # Request estimates using each force field, storing the request
    # object used to query the status of the results.
    for force_field_key in force_field_sources:

        force_field_source = force_field_sources[force_field_key]

        options = get_estimation_options(protocol_replacements.get(force_field_key, {}))

        requests[force_field_key] = estimator_client.request_estimate(property_set=data_set,
                                                                      force_field_source=force_field_source,
                                                                      options=options)

    # Wait for the results.
    should_run = True
    finished_force_fields = []

    while should_run:

        sleep(60)

        for force_field_key in force_field_sources:

            if force_field_key in finished_force_fields:
                continue

            results = requests[force_field_key].results(False)

            if isinstance(results, PropertyEstimatorResult) and len(results.queued_properties) > 0:
                continue

            logger.info(f'The server has completed {force_field_key}.')

            # Save the result to file.
            save_results(force_field_key, results)
            finished_force_fields.append(force_field_key)

        if len(finished_force_fields) == len(force_field_sources):
            should_run = False


if __name__ == "__main__":
    main()
