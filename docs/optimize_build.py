#!/usr/bin/env python3
"""
Post-build optimization script for acados documentation.
This script optimizes the built HTML, CSS, and JS files for better performance.
"""

import os
import re
import sys
from pathlib import Path
import gzip
import shutil

def add_defer_to_scripts(html_file):
    """Add defer attribute to script tags for better performance."""
    with open(html_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Add defer to FontAwesome script (it's not critical for initial render)
    content = re.sub(
        r'<script src="([^"]*fontawesome[^"]*)"([^>]*)>',
        r'<script defer src="\1"\2>',
        content
    )
    
    # Add defer to non-critical theme scripts
    # Keep documentation_options and doctools synchronous as they're needed early
    for script in ['bootstrap.js', 'pydata-sphinx-theme.js']:
        content = re.sub(
            f'<script src="([^"]*{script}[^"]*)"([^>]*)>',
            r'<script defer src="\1"\2>',
            content
        )
    
    # Add loading="lazy" to images
    # Exclude logo and important images from lazy loading
    content = re.sub(
        r'<img ([^>]*?)>',
        lambda m: f'<img {m.group(1)} loading="lazy">' if 'loading=' not in m.group(1) and 'logo' not in m.group(1).lower() and 'favicon' not in m.group(1).lower() else m.group(0),
        content
    )
    
    # Optimize external links - add rel="noopener noreferrer" for security and performance
    content = re.sub(
        r'<a ([^>]*?)href="https?://(?!docs\.acados\.org)([^"]*)"([^>]*?)>',
        lambda m: f'<a {m.group(1)}href="https://{m.group(2)}"{m.group(3)} rel="noopener noreferrer">' if 'rel=' not in m.group(1) and 'rel=' not in m.group(3) else m.group(0),
        content
    )
    
    with open(html_file, 'w', encoding='utf-8') as f:
        f.write(content)

def process_html_files(build_dir):
    """Process all HTML files in the build directory."""
    build_path = Path(build_dir)
    html_files = list(build_path.rglob('*.html'))
    
    print(f"Processing {len(html_files)} HTML files...")
    for html_file in html_files:
        add_defer_to_scripts(html_file)
    print("HTML optimization complete.")

def create_gzip_versions(build_dir):
    """Create gzip versions of static files for servers that support them."""
    build_path = Path(build_dir)
    extensions = ['.css', '.js', '.html', '.svg', '.json']
    
    files_to_compress = []
    for ext in extensions:
        files_to_compress.extend(build_path.rglob(f'*{ext}'))
    
    print(f"Creating gzip versions for {len(files_to_compress)} files...")
    compressed_count = 0
    for file_path in files_to_compress:
        # Skip already compressed files
        if file_path.suffix == '.gz':
            continue
        
        gz_path = Path(str(file_path) + '.gz')
        
        # Only compress if file is larger than 1KB
        if file_path.stat().st_size > 1024:
            with open(file_path, 'rb') as f_in:
                with gzip.open(gz_path, 'wb', compresslevel=9) as f_out:
                    shutil.copyfileobj(f_in, f_out)
            compressed_count += 1
    
    print(f"Created {compressed_count} gzip files.")

def main():
    # Determine build directory
    script_dir = Path(__file__).parent
    build_dir = script_dir / '_build'
    
    if not build_dir.exists():
        print(f"Build directory {build_dir} does not exist.")
        sys.exit(1)
    
    print("Starting post-build optimizations...")
    
    # Process HTML files to add defer/lazy loading
    process_html_files(build_dir)
    
    # Create gzip versions for better compression
    create_gzip_versions(build_dir)
    
    print("Post-build optimizations complete!")

if __name__ == '__main__':
    main()
