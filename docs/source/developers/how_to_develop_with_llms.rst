.. _developers-llm:

LLM-Assisted ImpactX Development
================================

Large Language Models (LLMs) can assist with ImpactX development tasks such as navigating the codebase, understanding existing implementations, writing new features, debugging, adding tests, and preparing pull requests.
This guide documents how the ImpactX repository is configured for LLM-based coding assistants and how to get the most out of them.

.. note::

   LLM-generated code should always be reviewed carefully.
   LLMs can hallucinate APIs, miss physics constraints, or produce code that compiles but is incorrect.
   The configurations described below help by providing accurate, up-to-date context, but developer judgment remains essential.

   The ImpactX community thus urges you to perform a careful, manual review of all LLM-generated code and documentation before asking for a review of your pull-request.
   This is important, otherwise you risk to waste the valuable time of our most proficient developers that will need to review your LLM-generated code.
   (Be considerate that ImpactX developers can prompt an LLM just as efficiently as you can. Your critical thinking skill to make sense of the LLM-generated code and make it sensible for review and maintainable for the long term is what is needed!)

.. note::

   This section is not understood as an endorsement of any of the listed (or unlisted) coding assistants or MCP services.
   Contributions to this section documenting further services, clients, skills, etc. are encouraged.

.. tip::

   When working with LLM coding assistants, keep in mind that *"most best practices are based on one constraint: [the] context window fills up fast, and performance degrades as it fills"* (`Claude Code Best Practices <https://code.claude.com/docs/en/best-practices>`__).
   Keep instructions concise, use plan modes, break complex tasks into focused sessions, and provide targeted context rather than overwhelming the assistant with information.


AGENTS.md / CLAUDE.md
---------------------

The repository can include an ``AGENTS.md`` file at its root (as well as a ``CLAUDE.md``, which directly points to ``AGENTS.md``).
These files are automatically loaded by LLM coding assistants (Claude Code reads ``CLAUDE.md``; other tools such as OpenAI Codex CLI read ``AGENTS.md``) to provide project-specific instructions.

The file contains in a compressed form instructions for an LLM agent:
With this file present, an LLM assistant working inside the ImpactX repository will automatically know how to build, test, and style code without being told each time.

To update these instructions, edit ``AGENTS.md``.
Keep this file under 300 lines to preserve LLM context.


Skills
------

LLM assistants such as Claude Code support reusable *skills* — scripted workflows that an assistant can execute on demand, automating multi-step tasks that follow a fixed procedure.
When defined for ImpactX, skills live in the ``.claude/skills/`` directory and use the prefix ``impactx-`` for easy discovery (start typing ``/impactx`` to see them).

To add a new skill, create a directory under ``.claude/skills/<skill-name>/`` containing a ``SKILL.md`` file that describes the step-by-step procedure.


Documentation Context via MCP Servers
--------------------------------------

LLM assistants work best when they can query up-to-date project documentation.
The :ref:`AI-Assisted Input File Design <ai_input_design>` workflow page describes how to set up `Model Context Protocol (MCP) <https://modelcontextprotocol.io>`__ servers for this purpose.
That setup is equally useful for development tasks: the same documentation context that helps write input files also helps an assistant understand ImpactX internals, AMReX, pyAMReX, and ABLASTR APIs, and conventions when writing C++ or Python code.

See :ref:`ai_input_design` for general MCP setup instructions with Context7.

ImpactX and Dependency Documentation on Context7
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`Context7 <https://context7.com>`__ is a service that indexes open-source project documentation and serves it through the MCP protocol.
Since ImpactX builds on top of `AMReX <https://amrex-codes.github.io/amrex/>`__, `pyAMReX <https://pyamrex.readthedocs.io>`__, `ABLASTR <https://github.com/BLAST-WarpX/warpx>`__ (shared with WarpX), and `openPMD-api <https://openpmd-api.readthedocs.io>`__, providing documentation for these dependencies alongside ImpactX documentation gives the assistant much richer context for development tasks.

The following documentation is available through Context7:

- **ImpactX**: `context7.com/blast-impactx/impactx <https://context7.com/blast-impactx/impactx>`__
- **WarpX / ABLASTR**: `context7.com/blast-warpx/warpx <https://context7.com/blast-warpx/warpx>`__
- **AMReX**: `context7.com/amrex-codes/amrex <https://context7.com/amrex-codes/amrex>`__
- **pyAMReX**: `context7.com/amrex-codes/pyamrex <https://context7.com/amrex-codes/pyamrex>`__
- **openPMD-api**: `context7.com/openpmd/openpmd-api <https://context7.com/openpmd/openpmd-api>`__
- **openPMD-viewer**: `context7.com/openpmd/openpmd-viewer <https://context7.com/openpmd/openpmd-viewer>`__
- **pybind11**: `context7.com/pybind/pybind11 <https://context7.com/pybind/pybind11>`__

When Context7 is connected, the assistant can look up any of those when needed:
AMReX data structures (e.g., ``MultiFab``, ``ParticleContainer``, ``Geometry``), pyAMReX and pybind11 binding patterns, ABLASTR utilities shared with WarpX, and openPMD I/O APIs directly, which is especially helpful when working, for instance, on:

- Beam-dynamics elements and lattice transport routines
- Space-charge and collective-effects solvers that use AMReX mesh data structures
- Particle routines built on ``amrex::ParticleContainer``
- Python bindings that wrap C++ classes via pybind11 and pyAMReX
- I/O and diagnostic code that interacts with AMReX plotfiles or openPMD

For instructions on configuring Context7 as an MCP server in your coding assistant (Claude Code, Cursor, VS Code, Windsurf, Codex CLI, and others), see the `Context7 client documentation <https://context7.com/docs/resources/all-clients>`__ and the :ref:`ai_input_design` page.
