from pathlib import Path
import re
from jinja2 import Environment


def github_path(toolname, scripts_tools):
    """Finds GitHub path for a tool in GRASS GIS with subtool handling"""

    auto_generated_indexes = [
        "database",
        "db",
        "display",
        "general",
        "imagery",
        "miscellaneous",
        "postscript",
        "raster",
        "raster3d",
        "temporal",
        "vector",
    ]

    # Special cases for documentation pages that don't match the pattern matching scheme
    special_docs = [
        {"name": "projectionintro", "path": "doc"},
        {"name": "grass_database", "path": "doc"},
        {"name": "databaseintro", "path": "db"},
    ]

    # Exit early if toolname is empty
    if not toolname:
        return None

    # Handle special cases
    if toolname in [x["name"] for x in special_docs]:
        return next((x["path"] for x in special_docs if x["name"] == toolname), None)

    # Convert filter() results to a list
    tool_matches = list(filter(lambda x: toolname in x, scripts_tools))

    # print(f"Tool Matches: {tool_matches}")

    # If there's exactly one match, return it immediately
    if len(tool_matches) == 1:
        return tool_matches[0]

    # Prioritize exact matches
    exact_match = next((x for x in tool_matches if x.endswith(toolname)), None)
    # print(f"Exact - Tool Matches: {exact_match}")
    if toolname in auto_generated_indexes:
        # print(f"Exact Match - Autogenerated: {exact_match}")
        return None
    if exact_match:
        return exact_match

    # Prefer major categories (raster, vector, imagery, temporal, etc.)
    category_match = next(
        (x for x in tool_matches if not x.startswith("scripts/")), None
    )
    # print(f"Category - Tool Matches: {category_match}")
    if category_match:
        return category_match

    # Check for subtools in the same directory (e.g., r.sim -> r.sim.water, r.sim.sediment)
    subtool_prefix = toolname + "."  # Ensure it matches the full prefix
    sub_tool_matches = [x for x in tool_matches if x.startswith(subtool_prefix)]

    # print(f"Subtool - Tool Matches: {sub_tool_matches}")
    # If exactly one subtool match exists, return it
    if len(sub_tool_matches) == 1:
        return sub_tool_matches[0]

    # If multiple subtools exist, prioritize based on name length (more specific subtools come first)
    if sub_tool_matches:
        sub_tool_matches.sort(key=len)  # Shorter names first
        return sub_tool_matches[0]  # Return the most specific match

    # Handle special case for "intro" pages
    # print(f"Intro Doc - Matches: {tool_matches}")
    if toolname.endswith("intro"):
        return toolname.replace("intro", "")

    # Special case for gui/wxpython/docs
    if toolname.startswith("g."):
        # print(f"GUI Docs: {toolname}")
        tool_dir = toolname.split(".")[-1]
        return f"gui/wxpython/{tool_dir}/{toolname}"

    # If nothing else, return the first available match
    return tool_matches[0] if tool_matches else None


def on_env(env: Environment, config, files):
    """Enable loopcontrols extension in Jinja2"""
    env.add_extension("jinja2.ext.loopcontrols")
    env.globals["github_path"] = github_path
    return env


def on_config(config):
    """
    Read the list of tools from the source directory and
    store it in MkDocs extra config. These are used to generate
    the correct links to the documentation in GitHub.
    """
    scripts_dir = Path("source")
    scripts_tools = []
    url_pattern = re.compile(
        r"https://github.com/OSGeo/grass/tree/main/([^ )]+)"
    )  # Read the mkdocs.yml file

    for file in scripts_dir.glob("*.md"):
        with file.open() as f:
            for line in f:
                match = url_pattern.search(line)
                if match:
                    toolname = match.group(1).strip()
                    scripts_tools.append(toolname)

    # Store in MkDocs extra config
    config["extra"]["scripts_tools"] = scripts_tools
    return config
