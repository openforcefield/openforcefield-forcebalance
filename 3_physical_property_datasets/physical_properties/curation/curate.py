import os

from nistdataselection.curatedataset import curate_data_set
from nistdataselection.generatereport import generate_report
from nistdataselection.parserawdata import parse_raw_data
from nistdataselection.utils.utils import SubstanceType

from propertyestimator.datasets import PhysicalPropertyDataSet
from propertyestimator.properties import ExcessMolarVolume, EnthalpyOfMixing, Density, DielectricConstant, \
    EnthalpyOfVaporization

# The set of those molecules which appeared in the VdW refit.
test_set_smiles = [
    'CC#N',
    'c1ccc(cc1)Cl',
    'c1ccsc1',
    'C(=O)O',
    'C(Cl)(Cl)Cl',
    'c1cc(cc(c1)Cl)Cl',
    'c1cc(cc(c1)Br)Br',
    'CCI',
    'C1CCNCC1',
    'CCOC(OCC)OCC',
    'c1ccc(cc1)F',
    'C(Br)Br',
    'C(Cl)Br',
    'CCOC(=O)CC(=O)C',
    'CC(=O)CCC(=O)O',
    'C(C(C(F)(F)F)(F)F)O',
    'COc1cccc(c1)Br',
    'C=CC#N',
    'CNCCCN',
    'CCn1ccnc1',
    'CCCI',
    'CCCCCCCCCCCCS',
    'COc1ccccc1O',
    'c1cc(cc(c1)I)F',
    'c1cc(ccc1F)I',
    'c1ccc2c(c1)ncs2',
    'Cc1ccc2c(c1)OCO2',
    'c1ccc(cc1)SCN=[N+]=[N-]',
    'c1cscc1C#N',
    'c1cc(sc1)C#N'
]

# The VdW parameters which were optimised.
optimized_vdw_smirks = [
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


def collate_nist_data(archive_directory, output_directory):
    """Collates a directory of NIST ThermoML archive files into
    more readily manipulable pandas csv files.
    """

    parse_raw_data(directory=archive_directory,
                   output_directory=output_directory,
                   retain_values=True,
                   conda_environment='benchmarking')


def curate(processed_data_directory, output_path):
    """Curate the benchmarking data set.
    """

    # Define the desired number of unique substances which should have data points
    # for each of the properties of interest
    desired_substances_per_property = {
        (EnthalpyOfMixing, SubstanceType.Binary): 10,
        (ExcessMolarVolume, SubstanceType.Binary): 10,
        (Density, SubstanceType.Pure): 30,
        (DielectricConstant, SubstanceType.Pure): 30,
        (EnthalpyOfVaporization, SubstanceType.Pure): 30
    }

    full_data_set = PhysicalPropertyDataSet()

    # Define the order of preference for which data binary substances should have.
    mixture_property_order = [
        [
            # We prioritise those molecules for which we have both binary enthalpies
            # of mixing and excess molar volumes.
            (ExcessMolarVolume, SubstanceType.Binary),
            (EnthalpyOfMixing, SubstanceType.Binary)
        ],
        [
            # Failing that, we pick molecules for which we only have enthalpies
            # of mixing.
            (EnthalpyOfMixing, SubstanceType.Binary)
        ],
        [
            # Finally, choose molecules for which we only have excess molar volumes.
            (ExcessMolarVolume, SubstanceType.Binary),
        ]
    ]

    # Build the mixture data sets.
    mixture_data_set = curate_data_set(processed_data_directory,
                                       mixture_property_order,
                                       desired_substances_per_property,
                                       required_smiles_to_include=None,
                                       smiles_to_exclude=[*test_set_smiles, 'O'],
                                       vdw_smirks_to_exercise=optimized_vdw_smirks,
                                       output_data_set_path='mixture_data_set.json')

    # We explicitly ask for aqueous mixture data.
    water_mixture_data_set = curate_data_set(processed_data_directory,
                                             mixture_property_order,
                                             desired_substances_per_property,
                                             required_smiles_to_include=['O'],
                                             smiles_to_exclude=test_set_smiles,
                                             vdw_smirks_to_exercise=optimized_vdw_smirks,
                                             output_data_set_path='water_mixture_data_set.json')

    full_data_set.merge(mixture_data_set)
    full_data_set.merge(water_mixture_data_set)

    # Next, build the pure data sets. Start by collating all of the previously chosen
    # molecules
    chosen_mixture_smiles = set()

    for properties in full_data_set.properties.values():

        for physical_property in properties:
            chosen_mixture_smiles.update([component.smiles for component in physical_property.substance.components])

    # Define the order of preference for which data pure substances should have.
    pure_property_order = [
        [
            (Density, SubstanceType.Pure),
            (DielectricConstant, SubstanceType.Pure)
        ],
        [
            (Density, SubstanceType.Pure)
        ],
        [
            (DielectricConstant, SubstanceType.Pure)
        ],
        [
            (EnthalpyOfVaporization, SubstanceType.Pure)
        ]
    ]

    # Ideally choose molecules for which we have also chosen binary data.
    # We exclude water as we did not aim to refit that in this release.
    pure_data_set = curate_data_set(processed_data_directory,
                                    pure_property_order,
                                    desired_substances_per_property,
                                    required_smiles_to_include=chosen_mixture_smiles,
                                    smiles_to_exclude=[*test_set_smiles, 'O'],
                                    vdw_smirks_to_exercise=optimized_vdw_smirks,
                                    minimum_data_points_per_property_per_smirks=3,
                                    output_data_set_path='pure_data_set_binary_compounds.json')

    chosen_pure_smiles = set()

    for properties in pure_data_set.properties.values():

        for physical_property in properties:
            chosen_pure_smiles.update([component.smiles for component in physical_property.substance.components])

    # Relax the criteria to include other molecules (again excluding water).
    pure_data_set.merge(curate_data_set(processed_data_directory,
                                        pure_property_order,
                                        desired_substances_per_property,
                                        required_smiles_to_include=None,
                                        smiles_to_exclude=[*test_set_smiles, 'O', *chosen_pure_smiles],
                                        vdw_smirks_to_exercise=optimized_vdw_smirks,
                                        minimum_data_points_per_property_per_smirks=3,
                                        output_data_set_path='pure_data_set.json'))

    full_data_set.merge(pure_data_set)

    with open(output_path, 'w') as file:
        file.write(full_data_set.json())


def report(data_set_path, report_name):
    """Generate a pdf report of the chosen data set.
    """

    generate_report(data_set_path, report_name, optimized_vdw_smirks)


def main():
    """An example method which calls each of the curation steps.
    """
    home_directory = os.path.expanduser("~")
    thermoml_archive_directory = os.path.join(home_directory, 'checked_thermoml_files')

    processed_data_directory = 'processed_data'
    collate_nist_data(thermoml_archive_directory, processed_data_directory)

    data_set_path = 'curated_data_set.json'
    curate(processed_data_directory, data_set_path)

    report(data_set_path, 'release_1_benchmark_set')


if __name__ == '__main__':
    main()
