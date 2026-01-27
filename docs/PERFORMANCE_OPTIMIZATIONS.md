# Documentation Performance Optimizations

This document describes the performance optimizations implemented for the acados documentation to improve the Google Lighthouse score.

## Optimizations Implemented

### 1. Custom CSS and JavaScript
- **File**: `_static/custom.css` - Optimizes font loading with `font-display: swap`
- **File**: `_static/custom.js` - Adds lazy loading for images
- These files are automatically included in the build via `conf.py`

### 2. Deferred Script Loading
- **File**: `optimize_build.py` - Post-build script that:
  - Adds `defer` attribute to non-critical JavaScript files (FontAwesome, Bootstrap, theme scripts)
  - Adds `loading="lazy"` to images for lazy loading
  - Creates gzip-compressed versions of static assets

### 3. Cache Control Headers
- **File**: `_headers` - Netlify-style headers file for cache control
- **File**: `_static/cache-config.txt` - Documentation for various web server cache configurations
- Static assets (CSS, JS, images, fonts): Cached for 1 year
- HTML files: Cached for 1 hour with must-revalidate

### 4. Build Process
The optimizations are automatically applied during the build process:
```bash
make html  # Builds HTML and runs optimize_build.py
```

## Performance Impact

These optimizations target the main Lighthouse performance issues:

1. **Reduce Unused CSS/JS**: 
   - Deferred loading of non-critical scripts
   - FontAwesome loaded with `defer` to prevent render blocking

2. **Efficient Cache Policy**:
   - Long cache times for static assets (1 year)
   - Shorter cache for HTML (1 hour) to ensure content freshness
   - Support for pre-compressed gzip files

3. **Lazy Loading**:
   - Images load only when needed
   - Reduces initial page load time

## Deployment Considerations

### For Netlify or Similar Platforms
The `_headers` file in the build output will automatically configure cache headers.

### For Apache
Use the directives from `_static/cache-config.txt` in your `.htaccess` or virtual host configuration.

### For Nginx
Add the location blocks from `_static/cache-config.txt` to your nginx configuration.

### For GitHub Pages
GitHub Pages automatically caches static assets, but you may want to configure additional optimizations through the repository settings.

## Testing

To test the optimizations locally:

1. Build the documentation:
   ```bash
   cd docs
   make clean
   make default  # or make doxygen && make html
   ```

2. Serve the built documentation:
   ```bash
   python -m http.server --directory _build 8000
   ```

3. Run Lighthouse audit on http://localhost:8000

## Maintenance

- The `optimize_build.py` script is run automatically during the build
- Custom CSS/JS files are small and focused on performance
- Cache headers should be reviewed periodically to ensure they align with deployment platform best practices
