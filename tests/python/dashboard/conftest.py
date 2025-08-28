import os
import time
import pytest
from seleniumbase import SB

from .utils import (
    DashboardTester,
    start_dashboard,
    wait_for_interaction_ready,
    wait_for_server_ready,
)


@pytest.hookimpl(hookwrapper=True, tryfirst=True)
def pytest_runtest_makereport(item, call):
    """Attach test phase reports to the node for later inspection.

    Allows fixtures to check failures in teardown via request.node.rep_call.
    """
    outcome = yield
    rep = outcome.get_result()
    setattr(item, "rep_" + rep.when, rep)


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
def reset_dashboard_inputs(dashboard, request):
    """
    Resets the dashboard to its default state before each test.

    This ensures all tests start from a clean, consistent baseline.
    """
    dashboard.sb.click("#Input_route")
    dashboard.sb.click("#reset_all_inputs_button")
    # Teardown: on failure, save a screenshot for debugging
    yield
    rep = getattr(request.node, "rep_call", None)
    if rep and rep.failed:
        try:
            os.makedirs("screenshots", exist_ok=True)
            ts = time.strftime("%Y%m%d-%H%M%S")
            name = f"{request.node.name}_{ts}.png"
            path = os.path.abspath(os.path.join("screenshots", name))
            # Use the underlying WebDriver to capture the screen
            dashboard.sb.driver.save_screenshot(path)
            print(f"Saved failed test screenshot: {path}")
        except Exception:
            pass
