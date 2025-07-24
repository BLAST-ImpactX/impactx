import pytest
from seleniumbase import SB

from .utils import (
    DashboardTester,
    start_dashboard,
    wait_for_interaction_ready,
    wait_for_server_ready,
)


@pytest.fixture(scope="session")
def dashboard():
    """
    Sets up a single ImpactX dashboard server and headless browser for all tests in the session.
    Automatically shuts down the server after all tests complete.
    """

    app_process = None

    try:
        app_process = start_dashboard()
        wait_for_server_ready(app_process)

        with SB(headless=True) as sb:
            sb.open("http://localhost:8080/index.html#/Input")
            wait_for_interaction_ready(sb)
            yield DashboardTester(sb)
    finally:
        if app_process:
            app_process.terminate()


@pytest.fixture(autouse=True)
def reset_dashboard_inputs(dashboard):
    """
    Resets the dashboard to its default state before each test.

    This ensures all tests start from a clean, consistent baseline.
    """
    dashboard.sb.click("#Input_route")
    dashboard.sb.click("#reset_all_inputs_button")
