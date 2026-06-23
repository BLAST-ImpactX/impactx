window.getSimStatusColor = function(status) {
    if (status === 'Completed') return 'success';
    if (status === 'In Progress') return 'info';
    if (status === 'Cancelled') return 'warning';
    if (status === 'Failed') return 'error';
    return 'default';
}

window.formatDate = function(isoString) {
    return new Date(isoString).toLocaleDateString('en-US', {
      month: 'short',
      day: 'numeric',
      year: 'numeric'
    });
  };

window.formatTime = function(isoString) {
return new Date(isoString).toLocaleTimeString([], {
    hour: '2-digit',
    minute: '2-digit',
    second: '2-digit'
});
};
