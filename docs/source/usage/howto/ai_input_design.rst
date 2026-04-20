.. _ai_input_design:

AI (LLM)-Assisted Input File Design
===================================

Large Language Models (LLMs) can accelerate the process of creating and modifying ImpactX Python scripts and input files.
By providing an LLM-based coding assistant with ImpactX documentation as context, users can describe an accelerator lattice and beam configuration in natural language and receive a draft Python script or parameter list file, ask for explanations of existing parameters, or request modifications to an existing configuration.

This workflow is equally applicable to:

- **ImpactX Python scripts** (:ref:`Parameters: Python <usage-python>`)
- **ImpactX parameter list files** (:ref:`Parameters: Inputs File <running-cpp-parameters>`)

.. note::

   LLM-generated Python scripts and input files should always be reviewed carefully before use.
   LLMs tend to hallucinate options that do not exist or might miss the mental context of your physics.
   The guide below tries to improve this situation by providing LLMs all ImpactX examples and documentation automatically, but this remains an active and rapidly evolving field of tooling.

   Validate physics parameters (lattice element lengths and strengths, beam energy and distribution, space-charge and CSR settings, integrator step counts) against your domain knowledge and the `ImpactX documentation <https://impactx.readthedocs.io>`__.


.. note::

   This section is not understood as an endorsement of any of the listed (or unlisted) coding assistants or MCP services.
   Contributions to this section documenting further services, clients, skills, etc. are encouraged.

How It Works: MCP Servers
-------------------------

`Model Context Protocol (MCP) <https://modelcontextprotocol.io>`__ servers are a standardized way to provide external context, such as library documentation, to LLM-based coding assistants.
When an MCP server is configured, the assistant can query up-to-date ImpactX documentation on demand, rather than relying solely on its training data.

Setting Up Context7 as an MCP Server
------------------------------------

`Context7 <https://context7.com>`__ is a service that indexes open-source project documentation and serves it through the MCP protocol.
ImpactX documentation is available at:

    `context7.com/blast-impactx/impactx <https://context7.com/blast-impactx/impactx>`__

When writing Python scripts, users will also encounter `pyAMReX <https://pyamrex.readthedocs.io>`__ APIs (e.g., for accessing mesh and particle data, or for extending simulations with callbacks):

- **pyAMReX**: `context7.com/amrex-codes/pyamrex <https://context7.com/amrex-codes/pyamrex>`__

To configure Context7 for your coding assistant, see the `Context7 documentation <https://context7.com/docs/resources/all-clients>`__.
Once connected, a coding assistant (Claude Code, Cursor, VS Code Copilot, Windsurf, etc.) can retrieve relevant sections of all the above projects in real time when helping you write or debug Python scripts and input files.

.. tip::

   Some clients support optional API key authentication.
   See the `Context7 client docs <https://context7.com/docs/resources/all-clients>`__ for details on API keys and OAuth options.

After creating a personal account, Context7 is free within limits for open source projects.

Best Practices
--------------

When working with LLM coding assistants, keep in mind that *"most best practices are based on one constraint: [the] context window fills up fast, and performance degrades as it fills"* (`Claude Code Best Practices <https://code.claude.com/docs/en/best-practices>`__).
Starting from examples and iterating incrementally as described below and making a plan with the LLM assistant helps keep sessions focused and productive.

#. **Start from examples.**
   Run your coding agent inside the ImpactX source directory.
   Point the assistant to an existing Python script or input file from the ``examples/`` directory (e.g., ``fodo``, ``chicane``, ``iota_lens``) and ask it to adapt the setup to your needs.

#. **Iterate incrementally.**
   Start with a minimal working lattice, verify it runs, then ask the assistant to add complexity (additional elements, space-charge or CSR, diagnostics, alignment errors, etc.).

#. **Cross-check with documentation.**
   Use the assistant to look up parameter meanings and valid ranges in the ImpactX docs, especially for numerical parameters like the integrator step counts (``nslice``), space-charge mesh resolution, and Poisson solver settings.

#. **Validate before production runs.**
   Always run a short test simulation and inspect the output before committing to a large-scale run.
   The ``examples/`` directory contains analysis scripts (``analysis_*.py``) that can serve as templates for validation.
