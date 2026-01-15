#!/usr/bin/env python3
#
# Copyright (c) The acados authors.
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#

"""
Generate sitemap.xml for acados documentation.

This script scans the built HTML documentation and creates a sitemap.xml file
that lists all pages with their last modification times and priority.
"""

import os
import datetime
from pathlib import Path
from xml.etree import ElementTree as ET
from xml.dom import minidom


# Files to exclude from the sitemap
SKIP_FILES = ['genindex.html', 'search.html', 'searchindex.html']


def get_last_modified(filepath):
    """Get the last modified time of a file."""
    timestamp = os.path.getmtime(filepath)
    return datetime.datetime.fromtimestamp(timestamp).strftime('%Y-%m-%d')


def prettify_xml(elem):
    """Return a pretty-printed XML string for the Element."""
    rough_string = ET.tostring(elem, encoding='utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="  ", encoding='utf-8').decode('utf-8')


def calculate_priority(url_path):
    """Calculate priority based on URL depth and importance."""
    # Root page has highest priority
    if url_path == 'index.html' or url_path == '':
        return '1.0'
    
    # Main section pages (one level deep)
    depth = url_path.count('/')
    if depth == 1:
        return '0.8'
    elif depth == 2:
        return '0.6'
    else:
        return '0.5'


def generate_sitemap(build_dir, base_url, output_file):
    """
    Generate sitemap.xml for the documentation.
    
    Args:
        build_dir: Path to the built HTML documentation directory
        base_url: Base URL of the documentation site (e.g., https://docs.acados.org)
        output_file: Path where sitemap.xml should be written
    """
    build_path = Path(build_dir)
    
    if not build_path.exists():
        print(f"Error: Build directory {build_dir} does not exist.")
        print("Please build the documentation first using 'make html'")
        return False
    
    # Create root element
    urlset = ET.Element('urlset')
    urlset.set('xmlns', 'http://www.sitemaps.org/schemas/sitemap/0.9')
    
    # Find all HTML files
    html_files = sorted(build_path.glob('**/*.html'))
    
    if not html_files:
        print(f"Warning: No HTML files found in {build_dir}")
        return False
    
    # Process each HTML file
    url_count = 0
    for html_file in html_files:
        # Get relative path from build directory
        rel_path = html_file.relative_to(build_path)
        
        # Skip certain files
        if html_file.name in SKIP_FILES:
            continue
        
        # Create URL entry
        url = ET.SubElement(urlset, 'url')
        
        # Add location
        loc = ET.SubElement(url, 'loc')
        url_path = str(rel_path).replace('\\', '/')
        loc.text = f"{base_url.rstrip('/')}/{url_path}"
        
        # Add last modified date
        lastmod = ET.SubElement(url, 'lastmod')
        lastmod.text = get_last_modified(html_file)
        
        # Add change frequency
        changefreq = ET.SubElement(url, 'changefreq')
        changefreq.text = 'weekly'
        
        # Add priority
        priority = ET.SubElement(url, 'priority')
        priority.text = calculate_priority(url_path)
        
        url_count += 1
    
    # Write to file
    xml_string = prettify_xml(urlset)
    
    # Remove the XML declaration line added by minidom and add our own
    lines = xml_string.split('\n')
    lines = [line for line in lines if not line.strip().startswith('<?xml')]
    xml_string = '<?xml version="1.0" encoding="UTF-8"?>\n' + '\n'.join(lines)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(xml_string)
    
    print(f"Successfully generated sitemap with {url_count} URLs")
    print(f"Sitemap written to: {output_file}")
    
    return True


def main():
    """Main entry point."""
    # Default paths
    script_dir = Path(__file__).parent
    build_dir = script_dir / '_build'
    output_file = build_dir / 'sitemap.xml'
    base_url = 'https://docs.acados.org'
    
    print("Generating sitemap.xml for acados documentation...")
    print(f"Build directory: {build_dir}")
    print(f"Base URL: {base_url}")
    print(f"Output file: {output_file}")
    print()
    
    success = generate_sitemap(build_dir, base_url, output_file)
    
    if success:
        print("\nSitemap generation completed successfully!")
        return 0
    else:
        print("\nSitemap generation failed!")
        return 1


if __name__ == '__main__':
    import sys
    sys.exit(main())
