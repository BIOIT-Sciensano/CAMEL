import pydantic

from camel.app.core.cameltestsuite import CamelTestSuite
from camel.app.core.errors import InvalidParameterError
from camel.app.core.parameter import Parameter
from camel.app.tools.dummytool import DummyTool


class TestTool(CamelTestSuite):
    """
    Tests for the tool class.
    """

    def test_params(self) -> None:
        """
        Tests adding and removing parameters to the dummy tool.
        :return: None
        """
        # Add a single parameter
        dummy_tool = DummyTool()
        dummy_tool.clear_parameters()
        dummy_tool.update_parameters(param_a='abc')
        self.assertEqual(len(dummy_tool.params), 1)
        self.assertIn('param_a', dummy_tool.params)
        cmd = dummy_tool.build_command()
        option_param_a = dummy_tool.params['param_a'].option
        self.assertIn(option_param_a, cmd.command)

        # Remove the same parameter
        dummy_tool.update_parameters(param_b='123')
        dummy_tool.update_parameters(param_a=False)
        self.assertEqual(len(dummy_tool.params), 1)
        self.assertIn('param_b', dummy_tool.params)
        cmd = dummy_tool.build_command()
        self.assertNotIn(option_param_a, cmd.command)
        self.assertIn(dummy_tool.params['param_b'].option, cmd.command)

        # Update parameter value
        dummy_tool.update_parameters(param_b='456')
        self.assertEqual(len(dummy_tool.params), 1)
        self.assertIn('param_b', dummy_tool.params)
        self.assertEqual(dummy_tool.params['param_b'].value, '456')
        cmd = dummy_tool.build_command()
        self.assertIn(f"{dummy_tool.params['param_b'].option} 456", cmd.command)

        # Add a flag
        dummy_tool.update_parameters(flag_a=True)
        self.assertEqual(len(dummy_tool.params), 2)
        self.assertIn('flag_a', dummy_tool.params)
        cmd = dummy_tool.build_command()
        self.assertIn(dummy_tool.params['flag_a'].option, cmd.command)

    def test_param_init(self) -> None:
        """
        Tests parameter creation.
        :return: None
        """
        # Create a valid parameter
        param = Parameter(**{'name': 'my_param', 'option': '--param', 'mandatory': False})
        self.assertEqual(param.name, 'my_param')
        self.assertEqual(param.option, '--param')

        # Create an invalid parameter (name missing)
        with self.assertRaises(pydantic.ValidationError):
            Parameter(**{'option': '--param', 'mandatory': False})

        # Create an invalid parameter (invalid field)
        with self.assertRaises(pydantic.ValidationError):
            Parameter(**{'name': 'param', 'invalid': 'abc'})

    def test_mandatory_parameter(self) -> None:
        """
        Tests if an error is raised when a mandatory parameter is not set.
        :return: None
        """
        dummy_tool = DummyTool()

        # Run without the mandatory parameter set -> should raise an error
        dummy_tool.clear_parameters()
        with self.assertRaises(InvalidParameterError):
            dummy_tool.run(self.running_dir)

        # run with the mandatory parameter set -> should not raise an error
        dummy_tool.update_parameters(mandatory='abc')
        dummy_tool.run(self.running_dir)

    def test_get_param_value(self) -> None:
        """
        Tests the get_param_value function.
        :return: None
        """
        dummy_tool = DummyTool()
        dummy_tool.clear_parameters()

        # Regular parameter
        dummy_tool.update_parameters(param_a='abc')
        self.assertEqual(dummy_tool.get_param_value('param_a'), 'abc')
        dummy_tool.update_parameters(param_a='def')
        self.assertEqual(dummy_tool.get_param_value('param_a'), 'def')

        # Flag parameter
        dummy_tool.update_parameters(flag_a=True)
        self.assertEqual(dummy_tool.get_param_value('flag_a'), True)
        dummy_tool.update_parameters(flag_a=False)
        self.assertEqual(dummy_tool.get_param_value('flag_a'), False)
