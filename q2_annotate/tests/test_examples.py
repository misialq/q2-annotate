from rachis.plugin.testing import TestPluginBase


class TestUsageExample(TestPluginBase):
    package = "q2_annotate.tests"

    def test_usage_examples(self):
        self.execute_examples()
