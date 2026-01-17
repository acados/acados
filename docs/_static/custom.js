/**
 * Custom JavaScript for acados documentation
 * Performance optimizations
 */

// Defer non-critical scripts
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
})();
