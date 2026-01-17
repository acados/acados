/**
 * Custom JavaScript for acados documentation
 * Performance optimizations
 */

(function() {
    'use strict';
    
    // Mark lazy-loaded images as loaded when they appear
    if ('IntersectionObserver' in window) {
        const imageObserver = new IntersectionObserver((entries, observer) => {
            entries.forEach(entry => {
                if (entry.isIntersecting) {
                    const img = entry.target;
                    img.classList.add('loaded');
                    observer.unobserve(img);
                }
            });
        });

        document.addEventListener('DOMContentLoaded', () => {
            document.querySelectorAll('img[loading="lazy"]').forEach(img => {
                imageObserver.observe(img);
            });
        });
    }
    
    // Reduce layout shifts by setting dimensions early
    document.addEventListener('DOMContentLoaded', () => {
        // Add aspect ratio to images without explicit dimensions
        document.querySelectorAll('img:not([width]):not([height])').forEach(img => {
            // Wait for image to load if not already loaded
            if (img.complete && img.naturalWidth && img.naturalHeight) {
                img.style.aspectRatio = `${img.naturalWidth} / ${img.naturalHeight}`;
            } else {
                img.addEventListener('load', function() {
                    if (this.naturalWidth && this.naturalHeight) {
                        this.style.aspectRatio = `${this.naturalWidth} / ${this.naturalHeight}`;
                    }
                });
            }
        });
    });
})();

