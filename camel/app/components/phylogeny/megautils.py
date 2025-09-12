from camel.app.tools.mega.mltreeconstruction import MLTreeConstruction
from camel.app.tools.mega.modelselection import ModelSelection


class MEGAUtils:
    """
    This class contains utility functions to work with the MEGA tools.
    """

    @staticmethod
    def update_model_selection_parameters(model_selection: ModelSelection, missing_data: str, branch_swap: str,
                                          site_cov_cutoff: int, threads: int) -> None:
        """
        Utility function to set the parameters for the ModelSelection tool.
        :param model_selection: Model selection tool instance
        :param missing_data: Missing data parameter value
        :param branch_swap: Branch swap parameter value
        :param site_cov_cutoff: Site coverage cutoff
        :param threads: Number of threads to use
        :return: None
        """
        if missing_data == 'complete_deletion':
            model_selection.update_parameters(missing_data_treatment='Complete deletion')
        elif missing_data == 'use_all_sites':
            model_selection.update_parameters(missing_data_treatment='Use all sites')
        elif missing_data == 'partial_deletion':
            model_selection.update_parameters(
                missing_data_treatment='Partial deletion', site_coverage_cutoff=site_cov_cutoff)
        else:
            raise ValueError(f"Invalid parameter value for missing data: '{missing_data}'")
        model_selection.update_parameters(branch_swap_filter=branch_swap.title().replace('_', ' '), threads=threads)

    @staticmethod
    def update_tree_building_parameters(tree_building: MLTreeConstruction, model: str, rates: str, bootstraps: int,
                                        missing_data: str, site_cov_cutoff: int, ml_method: str,
                                        branch_swap: str, threads: int) -> None:
        """
        Updates the parameters of the tree building tool.
        :param tree_building: Tree building tool instance
        :param bootstraps: Number of bootstrap replications
        :param rates: Rates among sites
        :param missing_data: Missing data treatment
        :param site_cov_cutoff: Site coverage cutoff (for partial deletion)
        :param model: Selection model
        :param ml_method: ML method (e.g. SPR3)
        :param branch_swap: Branch swap filter
        :param threads: Number of treads to use
        :return: None
        """
        tree_building.update_parameters(
            bootstrap_replications=bootstraps, test_of_phylogeny='Bootstrap method',
            model=model,
            branch_swap_filter=branch_swap.title().replace('_', ' '),
            heuristic_method=ml_method.upper(),
            threads=threads
        )

        # Add rates parameter
        if rates == 'G+I':
            tree_building.update_parameters(rates_among_sites='G+I')
            tree_building.update_parameters(gamma_categories='5')
        elif rates == 'G':
            tree_building.update_parameters(rates_among_sites='G')
            tree_building.update_parameters(gamma_categories='5')
        elif rates == 'I':
            tree_building.update_parameters(rates_among_sites='I')
        elif rates == 'U':
            tree_building.update_parameters(rates_among_sites='U')
        else:
            raise ValueError(f"Invalid rates parameter: {rates}")

        # Add missing data parameter
        if missing_data == 'complete_deletion':
            tree_building.update_parameters(missing_data_treatment='Complete deletion')
        elif missing_data == 'use_all_sites':
            tree_building.update_parameters(missing_data_treatment='Use all sites')
        elif missing_data == 'partial_deletion':
            tree_building.update_parameters(
                missing_data_treatment='Partial deletion', site_coverage_cutoff=site_cov_cutoff)
        else:
            raise ValueError(f"Invalid missing data parameter: {missing_data}")
